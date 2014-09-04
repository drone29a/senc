(ns moc.core
  (:require [clojure.core.matrix :as mx]
            [schema.core :as sc]
            [schema.macros :as sm]
            [moc.estimate.mle :as mle]
            [clojure.java.io :as io])
  (:use [clojure.tools.cli :only [cli]]
        [moc.community :only [community obj]]
        [munge.io.matrix-ctf :only [load-matrix save-matrix]]
        [munge.io.data-frame :only [save-data-frame]]
        [munge.io.core :only [load-ids]]
        [moc.schema :only [Vec Mat ProbVec]]
        [notion.metric :only [kl-cat]]
        [moc.util :only [proportional select-rows]])
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

(sm/defn select-comms :- {sc/Int Vec}
  "Returns a map of community index to community param vector, 
  given an object index and the matrix of all communities' params."
  [obj-membs :- Mat
   ])

(sm/defn e-step :- Mat
  "Calculate new object membership probabilities.
  The select-comms fn is used to select the communities
  of which an object is a member."
  [select-comms :- (sm/=> {sc/Int Vec} Mat sc/Int)
   objs :- [sc/Int]
   comm-params :- Mat])

(sm/defn m-step :- Mat
  "Calculate new community topic distribution parameters.
  The select-objs fn can be used to specify which objects should
  represent a specific community. E.g., a clique vs. all objects
  with non-zero membership probabilty."
  [select-objs :- (sm/=> {sc/Int Vec} Mat)
   comm-props :- Mat
   obj-membs :- Mat])

(sm/defn run
  "Run the algorithm. Main calls this, so can you."
  [obj-feats :- Mat
   obj-membs :- Mat
   comm-props :- Mat
   max-iter :- sc/Int
   comm-estimator :- (sm/=> ProbVec Mat Mat Mat sc/Int)]
  (let [num-comms (mx/row-count comm-props)]
    (loop [membs obj-membs
           props comm-props
           iter-count 0]
      (println "Iteration: " iter-count)
      (if (= max-iter iter-count)      
        {:obj-membs membs
         :comm-props props}
        (let [new-membs (e-step objs props)
              new-props (m-step props new-membs)
              ;; This was used before for "new-props"
              ;; (mx/matrix (map (partial comm-estimator obj-feats membs props)
              ;;                 (range num-comms)))
              ]
          (recur new-membs
                 new-props
                 (inc iter-count)))))))

(sm/defn estimate-props :- Mat
  [obj-feats :- Mat
   core-obj-idxs :- [[sc/Int]]]
  (mx/matrix (map (fn [idxs]
                    (proportional (reduce mx/add (select-rows obj-feats idxs))))
                  core-obj-idxs)))

(sm/defn normalize-log-prob :- ProbVec
  "Normalize log probabilities and return as plain probabilities.
  Probabilities < 1e-10 are dropped to zero."
  [log-probs :- Vec]
  (let [epsilon 1e-10
        threshold (- (Math/log epsilon) (Math/log (mx/row-count log-probs)))
        max-prob (mx/emax log-probs)]
    (proportional (mx/emap (comp #(if (> threshold %) 0 (Math/exp %))
                                 #(- % max-prob))
                           log-probs))))

(sm/defn feature-log-prob
  "Calculate probability of features given probability of each dimension.
  This is a multinomial PMF with the multinomial coefficient removed for efficiency.
  The probability output is unscaled!"
  [props :- ProbVec
   feat-counts :- Vec]
  (mx/dot feat-counts (mx/emap mle/safe-log props)))

(sm/defn estimate-memb :- Vec
  [obj-feats :- Vec
   comm-params :- Mat]
  ;; Calculate probability of obj for each community separately
  (mx/emap )
  ;; Normalize probabilities

  )

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
                    obj-membs]} (load-data (format "%s/%s_objs-affinity.csv" input-path i)
                                           (format "%s/%s_objs-feats.mtx" input-path i)
                                           (format "%s/%s_objs-membership.mtx" input-path i))
                    init-est-props (estimate-props obj-feats cores)
                    est-props (run obj-feats
                                obj-membs
                                init-est-props
                                max-iter
                                (partial mle/estimate-comm cores)
                                estimate-memb ; estimate object membership, takes obj, returns memb. uses get-obj-comms
                                )]
        (save-data-frame (format "%s/%s_mle-props.df" out-dir i) nil (mx/rows est-props))
        (save-data-frame (format "%s/%s_avg-props.df" out-dir i) nil (mx/rows init-est-props))))))



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
