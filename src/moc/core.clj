(ns moc.core
  (:require [clojure.core.matrix :as mx]
            [schema.core :as sc]
            [schema.macros :as sm]
            [moc.estimate.mle :as mle]
            [clojure.java.io :as io]
            [clojure.set :as set]
            [clojure.core.matrix.impl.pprint :as mpp :refer [pm]])
  (:use [clojure.tools.cli :only [cli]]
        [moc.community :only [community obj]]
        [munge.io.matrix-ctf :only [load-matrix save-matrix]]
        [munge.io.data-frame :only [save-data-frame]]
        [munge.io.core :only [load-ids]]
        [moc.schema :only [Vec Mat ProbVec BinVec BinMat]]
        [notion.metric :only [kl-cat]]
        [moc.util :only [proportional select-rows selected-rows binary-vector b-or]])
  (:gen-class))

(mx/set-current-implementation :vectorz)

(defn load-data
  [cores-path
   obj-feats-path
   obj-membs-path]
  (let [min-feat-count 1e-10]
    {:cores (vec (map (comp vec #(map dec %)) (load-ids cores-path)))
     :obj-feats (mx/add (load-matrix obj-feats-path) min-feat-count)
     :obj-membs (load-matrix obj-membs-path)}))

(sm/defn restricted-objs :- Mat
  "A convenience function for restricting object feature matrix to
  those objects belonging to the community's seed subgraph."
  [comm-objs :- BinMat
   obj-feats :- Mat]
  )

(sm/defn restricted-comms :- Mat
  "A convenience function for restricting the community params
  matrix according to those communities influential to a set of objects."
  [obj-comms :- BinMat
   comm-props :- Mat
   obj-idxs :- #{sc/Int}]
  (let [row-indicator (apply b-or (select-rows obj-comms (vec obj-idxs)))]
    (selected-rows comm-props row-indicator)))

(sm/defn estimate-props :- Mat
  [obj-feats :- Mat
   core-objs-vecs :- BinMat]
  (mx/matrix (map (fn [core-objs]
                    (proportional (reduce mx/add (selected-rows obj-feats core-objs))))
                  core-objs-vecs)))

(sm/defn estimate-memb :- Vec
  [obj-feats :- Vec
   comm-params :- Mat]
  ;; Calculate probability of object from weighted counts for all communities
  (->> (map proportional (mx/columns comm-params))
       (map (sm/fn [term-count :- sc/Num
                    weight-col :- Vec]
              (mx/emul term-count weight-col))
            obj-feats)
       (mx/sparse)
       (mx/transpose)
       (mx/rows)
       (map mx/esum)
       (mx/sparse)
       (proportional)))

;; TODO: Use the BinMat obj-groups and selected-rows instead of a select-obj-comms fn
;; or make m-step use a similar function. Check which is cleaner.
(sm/defn e-step :- Mat
  "Calculate new object membership probabilities.
  The select-comms fn is used to select the communities
  of which an object is a member."
  [select-obj-comms :- (sm/=> Mat Mat #{sc/Int}) ; This is restricted-comms with obj-comms partially applied
   obj-feats :- Mat
   comm-params :- Mat]
  ;; TODO: Make cleaner? This is a little goofy, building up a seq of comm-param matrices, one for each object by index.
  (let [obj-comm-params (map (comp (partial select-obj-comms comm-params) hash-set)
                             (range (mx/row-count obj-feats)))]
    (mx/matrix (map estimate-memb (mx/rows obj-feats) obj-comm-params))))

;; phi
(sm/defn subgraph-membership :- ProbVec
  "Calculate the community membership weights for the group of objects making up the subgraph."
  [objs-membs :- Mat
   obj-prop-num-obs :- ProbVec]
  (mx/mmul obj-prop-num-obs objs-membs))

;; Theta_c!=i
(sm/defn other-mixed-props :- Vec
  "The subgraph-memb is almost a ProbVec but the element corresponding the community being estimated
  will be zero.

  The result will be not quite a ProbVec either, as the community selected for estimation is not contributing."
  ;; TODO: we could return a proportional vector? But what if all zeros in the trivial case of objects that are
  ;; only members of a single community? Best to treat this result as non-proportional for consistency?
  [subgraph-memb :- Vec
   subset-comm-props :- Mat]
  (mx/mmul subgraph-memb
           subset-comm-props))

(sm/defn estimate-comm-props :- ProbVec
  [obj-groups :- BinMat
   obj-membs :- Mat
   obj-feats :- Mat
   comm-props :- Mat
   comm-idx :- sc/Int]
  (let [seed-objs (mx/get-column obj-groups comm-idx) ; objs part of the selected community seed subgraph

        ;; comms which objs from selected comm may be part of
        subset-comm-props (selected-rows comm-props
                                         (apply b-or (selected-rows obj-groups seed-objs))) 
        seed-objs-feats (reduce mx/add (selected-rows obj-feats seed-objs))

        ;; contribution proportion based on number of features per-object
        obj-prop-num-obs (proportional (reduce mx/add (mx/columns (selected-rows obj-feats seed-objs))))

        ;; TODO: shouldn't comm-mix-weight match what's calculated below in subgraph-membership call? (IT DOES!)
        ;; The issue is when we zero out the phi element for the selected community then phi is no longer normalized
        ;; That's why we then have to normalize the result in other-mixed-props.
        ;; Do we actually want it to be unnormalized? According to your equation I think so.
        ;; Think about this though... Maybe it doesn't matter?
        ;; Remember this is like mixing two distributions, one is the selected,
        ;; the other is everything else collapsed into one.
        
        ;; the "averaged" phi for the selected community's seed subgraph, this is the complete phi vector
        subgraph-memb (subgraph-membership (selected-rows obj-membs seed-objs) obj-prop-num-obs)

        ;; TODO: is it correct to mix communities according to per-object membership
        ;; but weighting it based on total number of object features?
        other-props (other-mixed-props (mx/mset subgraph-memb comm-idx 0) subset-comm-props)]

    (comment    (clojure.pprint/pprint subgraph-memb)
                (clojure.pprint/pprint seed-objs-feats)
                (clojure.pprint/pprint other-props)
                (clojure.pprint/pprint (mx/mget subgraph-memb comm-idx)))

    (mle/estimate-mixed-mult seed-objs-feats
                             other-props
                             (mx/mget subgraph-memb comm-idx))))

(sm/defn m-step :- Mat
  "Calculate new community topic distribution parameters.
  The select-objs fn can be used to specify which objects should
  represent a specific community. E.g., a clique vs. all objects
  with non-zero membership probabilty."
  [obj-groups :- BinMat
   obj-membs :- Mat
   obj-feats :- Mat
   comm-props :- Mat]
  (mx/matrix (map (partial estimate-comm-props obj-groups obj-membs obj-feats comm-props)
                  (range (mx/row-count comm-props)))))

(sm/defn run
  "Run the algorithm. Main calls this, so can you.
  Arguments:
  obj-feats - observed feature values
  obj-groups - the subgraph groups an object may be affiliated with, does not change
  obj-membs - initial estimate of object membership to communities
  comm-props - initial estimate of community topic parameters
  max-iter - maximum number of iterations"
  [obj-feats :- Mat
   obj-groups :- BinMat
   obj-membs :- Mat
   comm-props :- Mat
   max-iter :- sc/Int]
  (let [num-comms (mx/row-count comm-props)]
    (loop [membs obj-membs
           props comm-props
           iter-count 0]
      (println "Iteration: " iter-count)
      (println "Comm. properties:")
      (println (pm props))
      (if (= max-iter iter-count)      
        {:obj-membs membs
         :comm-props props}
        (let [new-membs (e-step (partial restricted-comms obj-groups) obj-feats props)
              new-props (m-step obj-groups membs obj-feats props)]
          (recur new-membs
                 new-props
                 (inc iter-count)))))))

;; TODO: Add schemas for matrices to ensure dimensions are correct. For example, the obj-feats matrix should have
;; the same number of columns as the comm-params matrix. I.e., we want a poor man's dependent types from predicate
;; schemas.

(defn -main
  [& args]
  (let [[opts args banner] (cli args
                                ["-h" "--help" "Print this help message and exit" :default false :flag true]
                                ["-o" "--out-dir" "Output directory" :default (System/getProperty "user.dir")]
                                ["-m" "--max-iter" "Max iterations" :default 10 :parse-fn #(Integer/parseInt %)])]
    (when (or (:help opts) (not= 1 (count args)))
      (println banner)
      (System/exit 0))
    
    (let [input-path (nth args 0)
          max-iter (:max-iter opts)
          out-dir (:out-dir opts)] 
      (let [{:keys [cores
                    obj-feats
                    obj-membs]} (load-data (format "%s/objs-affinity.csv" input-path)
                                           (format "%s/objs-feats.mtx" input-path)
                                           (format "%s/objs-membership.mtx" input-path))
                    initial-comm-props (estimate-props obj-feats cores)
                    initial-obj-membs (estimate-memb obj-feats initial-comm-props)
                    obj-groups nil
                    est-props (run
                                obj-feats
                                ;; TODO: Along with cores we have an obj-groups matrix which is based on
                                ;; object distane to all cores. Need to read/create obj-groups
                                obj-groups
                                initial-obj-membs
                                initial-comm-props
                                max-iter)]
        (save-data-frame (format "%s/mle-props.df" out-dir) nil (mx/rows est-props))
        (save-data-frame (format "%s/avg-props.df" out-dir) nil (mx/rows initial-comm-props))))))



(defn score
  [true-props est-props]
  (map kl-cat (mx/rows true-props) (mx/rows est-props)))

(sm/defn run-sim
  "Run the algorithm. Main calls this, so can you.
  Used for estimating parameters for simulated data with known true params."
  [score :- (sm/=> sc/Num Vec)
   obj-feats :- Mat
   obj-membs :- Mat
   comm-props :- Mat
   max-iter :- sc/Int
   comm-estimator :- (sm/=> Vec Mat Mat Mat sc/Int)]
  (let [num-comms (mx/row-count comm-props)]
    (loop [membs obj-membs
           props comm-props
           iter-count 0]
      (println (score props))
      (if (= max-iter iter-count)      
        props
        (recur obj-membs
               (mx/matrix (map (partial comm-estimator obj-feats membs props)
                               (range num-comms)))
               (inc iter-count))))))

(defn load-data-sim
  [cores-path
   obj-feats-path
   obj-membs-path
   comm-params-path]
  (let [min-feat-count 1e-10]
    {:cores (vec (map (comp vec #(map dec %)) (load-ids cores-path)))
     :obj-feats (mx/add (load-matrix obj-feats-path) min-feat-count)
     :obj-membs (load-matrix obj-membs-path)
     :true-props (load-matrix comm-params-path)}))

;; TODO: This main was built for experiments with synthetic data.
(defn -main-sim
  [& args]
  (let [[opts args banner] (cli args
                                ["-h" "--help" "Print this help message and exit" :default false :flag true]
                                ["-o" "--out-dir" "Output directory" :default (System/getProperty "user.dir")]
                                ["-m" "--max-iter" "Max iterations" :default 10 :parse-fn #(Integer/parseInt %)])]
    (when (or (:help opts) (not= 2 (count args)))
      (println banner)
      (System/exit 0))
    
    (let [input-path (nth args 0)
          num-probs (Integer/parseInt (nth args 1))
          max-iter (:max-iter opts)
          out-dir (:out-dir opts)] 
      (doseq [i (range 1 (inc num-probs))]
        (println (format "Running problem %s..." i))
        (let [{:keys [cores
                      obj-feats
                      obj-membs
                      true-props]} (load-data (format "%s/%s_objs-affinity.csv" input-path i)
                                              (format "%s/%s_objs-feats.mtx" input-path i)
                                              (format "%s/%s_objs-membership.mtx" input-path i)
                                              (format "%s/%s_comm-params.mtx" input-path i))
                      init-est-props (estimate-props obj-feats cores)
                      est-props (run (partial score true-props)
                                     obj-feats
                                     obj-membs
                                     init-est-props
                                     max-iter
                                     (partial mle/estimate-comm cores))]
          (save-data-frame (format "%s/%s_mle-props.df" out-dir i) nil (mx/rows est-props))
          (save-data-frame (format "%s/%s_avg-props.df" out-dir i) nil (mx/rows init-est-props)))))))


(comment

  (def result (sm/with-fn-validation
                (run (partial score
                              (:true-props data))
                  (:obj-feats data)
                  (:obj-membs data)
                  (estimate-props (:obj-feats data) (:cores data))
                  2
                  (partial moc.estimate.mle/estimate-comm (:cores data)))))

  )
