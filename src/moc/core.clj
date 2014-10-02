(ns moc.core
  (:require [clojure.core.matrix :as mx]
            [schema.core :as sc]
            [schema.macros :as sm]
            [moc.estimate.mle :as mle]
            [clojure.java.io :as io]
            [clojure.set :as set]
            [clojure.core.matrix.impl.pprint :as mpp :refer [pm]]
            [notion.metric :refer [kl-cat]]
            [clojure.pprint :refer [pprint]]
            [munge.core])
  (:use [clojure.tools.cli :only [cli]]
        [moc.community :only [community obj]]
        [munge.io.matrix-mm :only [load-matrix save-matrix]]
        [munge.io.data-frame :only [save-data-frame]]
        [munge.io.core :only [load-ids]]
        [moc.schema :only [Vec Mat ProbVec BinVec BinMat]]
        [moc.util :only [proportional select-rows selected-rows
                         binary-vector b-or safe-log log-multicat]])
  (:import [mikera.matrixx.impl SparseRowMatrix SparseColumnMatrix]
           [mikera.vectorz.impl SparseIndexedVector SparseHashedVector ZeroVector]
           [mikera.vectorz.util DoubleArrays])
  (:gen-class))

(set! *warn-on-reflection* true)
(mx/set-current-implementation :vectorz)

(defn print-matrix-info
  [msg m]
  (println msg
           (type m)
           (mx/sparse? m)
           (when (or (seq? m) (instance? java.util.List m))
             (format "seq contains: %s" (type (first m)))))
  m)

(sm/defn new-sparse-row-matrix :- Mat
  [sparse-vecs :- [Vec]]
  (SparseRowMatrix/create ^"[Lmikera.vectorz.AVector;" (into-array mikera.vectorz.AVector sparse-vecs)))

(sm/defn new-sparse-column-matrix :- Mat
  [sparse-vecs :- [Vec]]
  (SparseColumnMatrix/create ^"[Lmikera.vectorz.AVector;" (into-array mikera.vectorz.AVector sparse-vecs)))

(sm/defn new-sparse-indexed-vector :- Vec
  [v :- [sc/Num]]
  (let [data (double-array v)
        n (count data)
        nz-inds (DoubleArrays/nonZeroIndices data 0 n)
        nnz (count nz-inds)
        nz-data (double-array nnz)]
    (dotimes [i nnz]
      (aset nz-data i (aget data (aget nz-inds i))))
    (SparseIndexedVector/wrap n nz-inds nz-data)))

(sm/defn new-sparse-matrix :- Mat
  [rows :- [[sc/Num]]]
  (new-sparse-row-matrix (map new-sparse-indexed-vector rows)))

(sm/defn round-to-zero :- Mat
  "Returns a new matrix (need to recompute sparse index) with 
  values below threshold set to 0."
  [threshold :- sc/Num
   m :- (sc/either Mat Vec)]
  ;; TODO: best speed assumes row matrix, fixable through current core.matrix API?
  (let [vs (if (mx/vec? m) [m] (mx/rows m))]
    (SparseRowMatrix/create ^java.util.List (map (fn [v] (when (instance? SparseIndexedVector v)
                                                           (.roundToZero ^SparseIndexedVector v 0.001)))
                                                 vs))))

(sm/defn restricted-objs :- Mat
  "A convenience function for restricting object feature matrix to
  those objects belonging to the community's seed subgraph."
  [objs-groups :- BinMat
   objs-feats :- Mat
   comm-idx :- sc/Int]
  (->> (mx/get-column objs-groups comm-idx)
       (selected-rows objs-feats)))

(sm/defn restricted-comms :- Mat
  "A convenience function for restricting the community params
  matrix according to those communities influential to a set of objects."
  [objs-comms :- BinMat
   comms-props :- Mat
   obj-idxs :- #{sc/Int}]
  (let [row-indicator (apply b-or (select-rows objs-comms (vec obj-idxs)))]
    (selected-rows comms-props row-indicator)))

(sm/defn estimate-props :- Mat
  [objs-feat-vals :- Mat
   groups-objs :- BinMat]
  ;; Was using this instead of doseq:
  (let [ncols (mx/column-count objs-feat-vals)]
    (new-sparse-row-matrix (pmap (fn [group-objs]
                                   (let [accum (SparseIndexedVector/createLength ncols)]
                                     (doseq [row (mx/rows (selected-rows objs-feat-vals group-objs))]
                                       (when-not (instance? ZeroVector row)
                                         (mx/add! accum row)))
                                     (-> accum proportional)))
                                 groups-objs))))
(sm/defn estimate-memb :- Vec
  "Note: use a filtered comm-params matrix if you want to force zero membership 
  in some communities."
  [obj-feats :- Vec
   comms-params :- Mat]
  ;; Calculate probability of object from weighted counts for all communities

  ;; TODO: we don't need all columns, how can we only do calculations with non-zero columns?
  
  (let [v (SparseIndexedVector/createLength (mx/row-count comms-params))
        ;;weight (sm/fn )
        

        ;; Was using this for testing the chubby cons
        ;; result (->> (mx/columns comms-params)
                    
        ;;             (print-matrix-info "after first prop")
                    
        ;;             (new-sparse-column-matrix))
        ]
    (->> (mx/columns comms-params)
         (map proportional)
         (map (sm/fn [term-count :- sc/Num
                      weight-col :- Vec]
                (mx/emul! weight-col term-count)
                weight-col)
              obj-feats)
         ;; TODO: can remove vec? removing vec makes it really slooow. coll -> array faster than seq -> array?
         (vec)
         (new-sparse-column-matrix)
         (mx/rows)
         ;; TODO: is there a map we can use that will generate a sparse vector?
         (map mx/esum)
         ;; TODO: can remove vec?
         (vec)
         (new-sparse-indexed-vector)
         (proportional))))

;; TODO: convert into a function that can be used for computing single object memberships?
(sm/defn estimate-membs :- Mat
  "Estimate the membership of a group of objects, specified by index."
  [objs-feat-vals :- Mat
   objs-groups :- BinMat
   comms-props :- Mat
   obj-idxs :- [sc/Int]]
  (let [num-objs (mx/row-count objs-feat-vals)
        num-comms (mx/row-count comms-props)
        membs (SparseRowMatrix/create num-objs num-comms)]
    ;; TODO: Fix this below?
    ;; Memory still grows (because of mapping a closure?) but
    ;; doesn't appear to OOM. Fills up heap to max and then does a big cleanup.
    ;; Hypothesis is that now it clean up the mess since we're not holding the
    ;; head. Still annoyed that it gets so big in first place....
    ;; An attempt to not hold the head...
    (doseq [[idx obj-row] (map-indexed vector (pmap (fn [obj-idx]
                                                  (estimate-memb (mx/get-row objs-feat-vals obj-idx)
                                                                 (restricted-comms objs-groups
                                                                                   comms-props
                                                                                   (hash-set obj-idx))))
                                                    obj-idxs))]
      (mx/set-row! membs idx obj-row))

    ;; Avoids the memory explosion..
    ;; (dotimes [obj-idx num-objs]
    ;;   (let [obj-row (estimate-memb (mx/get-row objs-feat-vals obj-idx)
    ;;                                 (restricted-comms objs-groups
    ;;                                                   comms-props
    ;;                                                   (hash-set obj-idx)))]
    ;;        ;; TODO: get this non-zeros and row setter under control. core.matrix should support this!!!
    ;;     (mx/set-row! membs obj-idx obj-row)))
    
    membs))

;; TODO: Use the BinMat objs-groups and selected-rows instead of a select-obj-comms fn
;; or make m-step use a similar function. Check which is cleaner.
(sm/defn e-step :- Mat
  "Calculate new object membership probabilities.
  The select-comms fn is used to select the communities
  of which an object is a member."
  [select-obj-comms :- (sm/=> Mat Mat #{sc/Int}) ; This is restricted-comms with obj-comms partially applied
   objs-feats :- Mat
   comms-params :- Mat]
  ;; TODO: Make cleaner? This is a little goofy, building up a seq of comm-param matrices, one for each object by index.
  ;; TODO: Does this map get executed in the consuming thread run from pmap or by main thread before going to pmap?
  (let [obj-comm-params (map (comp (partial select-obj-comms comms-params) hash-set)
                             (range (mx/row-count objs-feats)))]
    (new-sparse-row-matrix (pmap estimate-memb (mx/rows objs-feats) obj-comm-params))))

;; phi
(sm/defn subgraph-membership :- ProbVec
  "Calculate the community membership weights for the group of objects making up the subgraph."
  [objs-membs :- Mat
   obj-prop-num-obs :- ProbVec]
  (mx/mmul obj-prop-num-obs objs-membs))

;; Theta_c!=i
(sm/defn mixed-props :- Vec
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
  [groups-objs :- BinMat
   objs-groups :- BinMat
   objs-membs :- Mat
   objs-feats :- Mat
   comms-props :- Mat
   comm-idx :- sc/Int]
  (let [seed-objs (mx/get-row groups-objs comm-idx) ; objs part of the selected community seed subgraph

        ;; comms which objs from selected comm may be part of
        subset-comm-props (selected-rows comms-props
                                         (apply b-or (selected-rows objs-groups seed-objs)))

        ;; TODO: Need to update matrix_api.clj to use SparseIndexedVector.
        ;;       Or it was just that addCopy isn't supported but regular add is fine?
        seed-objs-feats (reduce mx/add! (SparseIndexedVector/createLength (mx/column-count objs-feats))
                                (selected-rows objs-feats seed-objs))

        ;; contribution proportion based on number of features per-object
        obj-prop-num-obs (proportional (reduce mx/add! (SparseIndexedVector/createLength (mx/row-count objs-feats))
                                               (mx/columns (selected-rows objs-feats seed-objs))))

        ;; TODO: shouldn't comm-mix-weight match what's calculated below in subgraph-membership call? (IT DOES!)
        ;; The issue is when we zero out the phi element for the selected community then phi is no longer normalized
        ;; That's why we then have to normalize the result in other-mixed-props.
        ;; Do we actually want it to be unnormalized? According to your equation I think so.
        ;; Think about this though... Maybe it doesn't matter?
        ;; Remember this is like mixing two distributions, one is the selected,
        ;; the other is everything else collapsed into one.
        
        ;; the "averaged" phi for the selected community's seed subgraph, this is the complete phi vector
        subgraph-memb (subgraph-membership (selected-rows objs-membs seed-objs) obj-prop-num-obs)

        ;; TODO: is it correct to mix communities according to per-object membership
        ;; but weighting is based on total number of object features?
        other-props (mixed-props (mx/mset subgraph-memb comm-idx 0) subset-comm-props)]
    (when (= comm-idx 2)      
      ;;(println (pm other-props)) ; this was all 0s
      ;;(println (pm subgraph-memb)) ; has some nonzeros
      ;; (println (pm (mx/mget subgraph-memb comm-idx))) ; it is 0. what?!??
      ;; (println (pm (mx/get-column objs-membs comm-idx))) ; we need to lookup specific cells
      ;; (println (pm obj-prop-num-obs)) ; has some nonzeros, need to compare cells
      (println (seq (mx/non-zero-indices seed-objs))))
    (mle/estimate-mixed-mult seed-objs-feats
                             other-props
                             (mx/mget subgraph-memb comm-idx))))

(sm/defn score-params :- Vec
  [obj-groups :- BinMat
   obj-feats :- Mat
   comm-props :- Mat]
  (new-sparse-indexed-vector (pmap log-multicat
                                   (mx/rows comm-props)
                                   (pmap (comp (partial reduce mx/add)
                                               (partial restricted-objs obj-groups obj-feats))
                                         (range (mx/row-count comm-props))))))

;; TODO: why did I name this with alt- prefix
;; TODO: pass groups-objs from -main -> run -> m-step -> here
;;       avoids having to call getColumn, which is slow.
(sm/defn alt-mixed-prop :- ProbVec
  "Create a mixture distribution given a set of objects, their membership weights,
  and community feature distributions."
  [objs-groups :- BinMat
   objs-feats :- Mat
   objs-membs :- Mat
   comms-props :- Mat
   comm-idx :- sc/Int]
  (let [seed-objs (mx/get-column objs-groups comm-idx) ; objs part of the selected community seed subgraph
        subset-comms-props (selected-rows comms-props
                                          (apply b-or (selected-rows objs-groups seed-objs)))
        ;; TODO: loop or doseq faster?
        seed-objs-feats (reduce mx/add!
                                (SparseIndexedVector/createLength (mx/column-count objs-feats))
                                (selected-rows objs-feats seed-objs))
        
        ;; contribution proportion based on number of features per-object
        obj-prop-num-obs (proportional (reduce mx/add!
                                               (SparseIndexedVector/createLength (mx/row-count objs-feats))
                                               (mx/columns (selected-rows objs-feats
                                                                          seed-objs))))
        ;; the complete phi vector
        subgraph-memb (subgraph-membership (selected-rows objs-membs seed-objs)
                                           obj-prop-num-obs)]

    ;; These are SparseIndexedVector and SparseRowMatrix
    (mx/mmul subgraph-memb subset-comms-props)))

(sm/defn m-step :- {:comms-props Mat
                    :changed sc/Bool
                    :change-count sc/Int}
  "Calculate new community topic distribution parameters.
  The select-objs fn can be used to specify which objects should
  represent a specific community. E.g., a clique vs. all objects
  with non-zero membership probabilty."
  [groups-objs :- BinMat
   objs-groups-possible :- BinMat
   objs-membs :- Mat
   objs-feats :- Mat
   comms-props :- Mat]
  ;; Keep only the improvements. Is this nonsensical?
  (let [scorer (partial score-params objs-groups-possible objs-feats)
        new-props (->> (new-sparse-row-matrix (pmap (partial estimate-comm-props
                                                             groups-objs
                                                             objs-groups-possible
                                                             objs-membs
                                                             objs-feats
                                                             comms-props)
                                                    (range (mx/row-count comms-props))))
                       (round-to-zero 1e-4))
        old-mixed-props (new-sparse-row-matrix (pmap (partial alt-mixed-prop
                                                              objs-groups-possible
                                                              objs-feats
                                                              objs-membs
                                                              comms-props)
                                                     (range (mx/row-count comms-props))))
        new-mixed-props (new-sparse-row-matrix (pmap (fn [comm-idx]
                                                       (alt-mixed-prop objs-groups-possible
                                                                       objs-feats
                                                                       objs-membs
                                                                       new-props
                                                                       comm-idx))
                                                     (range (mx/row-count comms-props))))
        changed (atom false)
        change-count (atom 0)]
    {:comms-props (new-sparse-row-matrix (pmap (fn [[old-score old-v] [new-score new-v]]
                                                 (if (>= old-score new-score)
                                                   old-v
                                                   (do (swap! changed not)
                                                       (swap! change-count inc)
                                                       new-v)))
                                               (map vector (scorer old-mixed-props) comms-props)
                                               (map vector (scorer new-mixed-props) new-props)))
     :changed @changed
     :change-count @change-count}))

(sm/defn run
  "Run the algorithm. Main calls this, so can you.
  Arguments:
  objs-feat-vals - observed feature values
  groups-objs - the seed objects of each group
  objs-groups - the subgraph groups an object may be affiliated with, does not change
  objs-membs - initial estimate of object membership to communities
  comms-props - initial estimate of community topic parameters
  max-iter - maximum number of iterations"
  [objs-feat-vals :- Mat
   groups-objs :- BinMat
   objs-groups :- BinMat
   objs-membs :- Mat
   comms-props :- Mat
   max-iter :- sc/Int]
  (let [num-comms (mx/row-count comms-props)]
    (loop [membs objs-membs
           props comms-props
           iter-count 0]

      (println (format "\nIteration: %d, at: %s" iter-count (str (java.util.Date.))))
      ;;(pprint (score-params objs-groups objs-feat-vals props))
      
      (if (= max-iter iter-count)
        (do (println "Maximum number of iterations reached.")
            {:objs-membs membs
             :comms-props props})
        (let [new-membs (e-step (partial restricted-comms objs-groups) objs-feat-vals props)
              {changed :changed
               change-count :change-count
               new-props :comms-props} (m-step groups-objs objs-groups membs objs-feat-vals props)]
          (if (not changed)
            (do (println "Community distributions have converged.")
                {:objs-membs membs
                 :comms-props props})
            (do (println (format "Number of communities changed: %d" change-count))
                (recur new-membs
                       new-props
                       (inc iter-count)))))))))

;; TODO: Add schemas for matrices to ensure dimensions are correct. For example, the objs-feat-vals matrix should have
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

    (println "Starting up...")
    
    (let [input-path (nth args 0)
          max-iter (:max-iter opts)
          out-dir (:out-dir opts)
          objs-feat-vals (load-matrix (format "%s/objs-feat-vals.mm" input-path))
          objs-groups (load-matrix (format "%s/objs-groups.mm" input-path))
          groups-objs (load-matrix (format "%s/groups-objs.mm" input-path))
          _ (println "Loaded input.")
          initial-comms-props (->> (estimate-props objs-feat-vals groups-objs)
                                   (round-to-zero 1e-4))
          _ (println "Initial community proportions found.")
          seed-objs-idx (->> (mx/non-zero-indices groups-objs)
                             (mapcat identity)
                             (distinct))
          initial-objs-membs (estimate-membs objs-feat-vals objs-groups initial-comms-props seed-objs-idx)
          _ (println "Initial object memberships found.")
          {est-comms-props :comms-props
           est-objs-membs :objs-membs} (run
                                         objs-feat-vals
                                         groups-objs
                                         objs-groups
                                         initial-objs-membs
                                         initial-comms-props
                                         max-iter)
          _ (println "Estimating all objects membership...")
          num-objs (mx/row-count objs-feat-vals)
          all-est-objs-membs (estimate-membs objs-feat-vals objs-groups est-comms-props (range num-objs))]
      
      (println "Saving output...")
      (save-matrix est-comms-props (format "%s/mle-comms-props.df" out-dir))
      (save-matrix all-est-objs-membs (format "%s/mle-objs-membs.df" out-dir)))))



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
                      true-props]} (load-data-sim (format "%s/%s_objs-affinity.csv" input-path i)
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
