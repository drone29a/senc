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

(defn write-log
  [msg]
  (println (format "%s - %s" (str (java.util.Date.)) msg)))

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

(sm/defn create-sparse-indexed-vector :- Vec
  [length :- sc/Int]
  (SparseIndexedVector/createLength length))

(sm/defn new-sparse-matrix :- Mat
  [rows :- [[sc/Num]]]
  (new-sparse-row-matrix (map new-sparse-indexed-vector rows)))

(sm/defn round-to-zero! :- Mat
  "Returns a new matrix (need to recompute sparse index) with 
  values below threshold set to 0."
  [threshold :- sc/Num
   m :- (sc/either Mat Vec)]
  ;; TODO: best speed assumes row matrix, fixable through current core.matrix API?
  (let [vs (if (mx/vec? m) [m] (mx/rows m))]
    (SparseRowMatrix/create ^java.util.List (map (fn [v] (when (instance? SparseIndexedVector v)
                                                           (-> (.roundToZero ^SparseIndexedVector v threshold)
                                                               (proportional))))
                                                 vs))))

(sm/defn restricted-objs :- Mat
  "A convenience function for restricting object feature matrix to
  those objects belonging to the community's seed subgraph."
  [seed-groups-objs :- BinMat
   objs-feats :- Mat
   comm-idx :- sc/Int]
  (->> (mx/get-row seed-groups-objs comm-idx)
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
  "Used to calculate the initial estimates for community feature proportions."
  [objs-feat-vals :- Mat
   groups-objs :- BinMat]
  ;; Was using this instead of doseq:
  (let [ncols (mx/column-count objs-feat-vals)
        result (new-sparse-row-matrix (pmap (fn [group-objs]
                                              (let [accum (SparseIndexedVector/createLength ncols)]
                                                (doseq [row (mx/rows (selected-rows objs-feat-vals group-objs))]
                                                  (when-not (instance? ZeroVector row)
                                                    (mx/add! accum row)))
                                                (-> accum proportional)))
                                            (mx/rows groups-objs)))]
    result))

;; TODO: pass in a column-matrix version of comm-sparams for speed up?
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
  "Estimate the membership of a group of objects, specified by index. Used to limit estimation
  to only those objects belonging to seed groups."
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
    (doseq [[idx obj-row] (pmap (fn [obj-idx]
                                  [obj-idx (estimate-memb (mx/get-row objs-feat-vals obj-idx)
                                                          (restricted-comms objs-groups
                                                                            comms-props
                                                                            (hash-set obj-idx)))])
                                obj-idxs)]
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
  (let [num-objs (mx/row-count objs-feats)
        obj-comm-params (map (comp (partial select-obj-comms comms-params) hash-set)
                             (range num-objs))
        ;; TODO: you were already estimating all object memberships?!
        result (new-sparse-row-matrix (pmap estimate-memb (mx/rows objs-feats) obj-comm-params))]
    result))

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
  (let [;;seed-objs (mx/get-row groups-objs comm-idx) ; objs part of the selected community seed subgraph
        num-objs (mx/column-count groups-objs)
        seed-objs (let [objs (mx/get-column objs-membs comm-idx)
                        objs-bin (create-sparse-indexed-vector num-objs)]
                    (doseq [nzi (mx/non-zero-indices objs)]
                      (mx/mset! objs-bin nzi 1))
                    objs-bin)

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

    (comment (when (= comm-idx 0)
               (println "estimate-comm-props")
               (println (mx/esum other-props)) ; this was all 0s ; almost 0
               (println (mx/esum subgraph-memb)) ; has some nonzeros ; sums to 0.69
               (println (pm (mx/mget subgraph-memb comm-idx))) ; it is 0. what?!?? it is 0 again, what?!?!
               (println (mx/esum (mx/get-column objs-membs comm-idx))) ; we need to lookup specific cells
               (println (mx/esum obj-prop-num-obs)) ; has some nonzeros, need to compare cells
               (println (seq (mx/non-zero-indices seed-objs)))))

    (mle/estimate-mixed-mult seed-objs-feats
                             other-props
                             (mx/mget subgraph-memb comm-idx))))

(sm/defn score-params :- Vec
  [groups-objs :- BinMat
   objs-feats :- Mat
   comms-props :- Mat]
  (new-sparse-indexed-vector (pmap log-multicat
                                   (mx/rows comms-props)
                                   (pmap (comp (partial reduce mx/add)
                                               (partial restricted-objs groups-objs objs-feats))
                                         (range (mx/row-count comms-props))))))

;; TODO: why did I name this with alt- prefix, should this and estimate-comm-props be merged???
;; It's everything but the selected?  alt -> alternative?
(sm/defn alt-mixed-prop :- ProbVec
  "Create a mixture distribution given a set of objects, their membership weights,
  and community feature distributions."
  [groups-objs :- BinMat
   objs-groups :- BinMat   
   objs-feats :- Mat
   objs-membs :- Mat
   comms-props :- Mat
   comm-idx :- sc/Int]
  (let [seed-objs (mx/get-row groups-objs comm-idx) ; objs part of the selected community seed subgraph
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

    (comment (when (= comm-idx 0)
               (println "alt-mixed-prop")
               (println (mx/esum subgraph-memb)) ; has some nonzeros ; sums to 0.69
               (println (pm (mx/mget subgraph-memb comm-idx))) ; it is 0. what?!?? it is 0 again, what?!?!
               (println (mx/esum (mx/get-column objs-membs comm-idx))) ; we need to lookup specific cells
               (println (mx/esum obj-prop-num-obs)) ; has some nonzeros, need to compare cells
               (println (seq (mx/non-zero-indices seed-objs)))))
    
    ;; These are SparseIndexedVector and SparseRowMatrix
    (mx/mmul subgraph-memb subset-comms-props)))

(sm/defn m-step :- {:comms-props Mat
                    :changed sc/Bool
                    :change-count sc/Int
                    :changed-ids #{sc/Int}}
  "Calculate new community topic distribution parameters.
  The groups-objs matrix can be used to specify which objects should
  represent a specific community. E.g., a clique vs. all objects
  with non-zero membership probabilty."
  [groups-objs :- BinMat
   objs-groups-possible :- BinMat
   objs-membs :- Mat
   objs-feats :- Mat
   comms-props :- Mat]
  ;; Keep only the improvements. Is this nonsensical?
  ;; We want to score using all objects since we estimated comms with only seed group objects.
  (let [scorer (partial score-params (new-sparse-row-matrix (mx/columns objs-groups-possible)) objs-feats)
        num-comms (mx/row-count comms-props)
        new-props (->> (new-sparse-row-matrix (pmap (partial estimate-comm-props
                                                             groups-objs
                                                             objs-groups-possible
                                                             objs-membs
                                                             objs-feats
                                                             comms-props)
                                                    (range num-comms)))
                       (round-to-zero! 1e-4))
        old-mixed-props (new-sparse-row-matrix (pmap (partial alt-mixed-prop
                                                              groups-objs
                                                              objs-groups-possible
                                                              objs-feats
                                                              objs-membs
                                                              comms-props)
                                                     (range num-comms)))
        new-mixed-props (new-sparse-row-matrix (pmap (partial alt-mixed-prop
                                                              groups-objs
                                                              objs-groups-possible
                                                              objs-feats
                                                              objs-membs
                                                              new-props)
                                                     (range num-comms)))
        changed (atom false)
        change-count (atom 0)
        changed-ids (atom #{})]
    {:comms-props (new-sparse-row-matrix (pmap (fn [[old-score old-v] [new-score new-v] comm-idx]
                                                 (if (>= old-score new-score)
                                                   old-v
                                                   (do (swap! changed not)
                                                       (swap! change-count inc)
                                                       (swap! changed-ids conj comm-idx)
                                                       new-v)))
                                               (map vector (scorer old-mixed-props) comms-props)
                                               (map vector (scorer new-mixed-props) new-props)
                                               (range num-comms)))
     :changed @changed
     :change-count @change-count
     :changed-ids @changed-ids}))

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
           iter-count 0
           last-changed-comms #{}
           all-changed-comms #{}]

      (write-log (format "Iteration: %d" iter-count))
      ;;(pprint (score-params groups-objs objs-feat-vals props))
      
      (if (= max-iter iter-count)
        (do (write-log "Maximum number of iterations reached.")
            {:objs-membs membs
             :comms-props props})
        (let [new-membs (e-step (partial restricted-comms objs-groups) objs-feat-vals props)
              {changed :changed
               change-count :change-count
               changed-ids :changed-ids
               new-props :comms-props} (m-step groups-objs objs-groups new-membs objs-feat-vals props)
              ;; TODO: make sure to check for NOTE s too...
              ;; NOTE: We were using past iteration membs before, but is that wrong?
              ;;       Converges faster when using current iteration's membs
              ;;(m-step groups-objs objs-groups membs objs-feat-vals props)
              ]
          (if (not changed)
            (do (write-log "Community distributions have converged.")
                {:objs-membs membs
                 :comms-props props})
            (do (write-log (format "Number of communities changed: %d" change-count))
                (write-log (format "Number of different communities from last iteration: %d" (count (set/difference changed-ids last-changed-comms))))
                (write-log (format "Number of first-time changed communities: %d" (count (set/difference changed-ids all-changed-comms))))
                (recur new-membs
                       new-props
                       (inc iter-count)
                       changed-ids
                       (set/union all-changed-comms changed-ids)))))))))

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

    (write-log "Starting up...")
    
    (let [input-path (nth args 0)
          max-iter (:max-iter opts)
          out-dir (:out-dir opts)
          objs-feat-vals (load-matrix (format "%s/objs-feat-vals.mm" input-path))
          objs-groups (load-matrix (format "%s/objs-groups.mm" input-path))
          groups-objs (load-matrix (format "%s/groups-objs.mm" input-path))
          _ (write-log "Loaded input.")
          initial-comms-props (->> (estimate-props objs-feat-vals groups-objs)
                                   (round-to-zero! 1e-4))
          _ (write-log "Initial community proportions found.")
          ;; All the objects in at least one seed group
          seed-objs-idx (->> (mx/non-zero-indices groups-objs)
                             (mapcat identity)
                             (distinct))
          ;;initial-objs-membs (estimate-membs objs-feat-vals objs-groups initial-comms-props seed-objs-idx)
          num-objs (mx/column-count groups-objs)
          initial-objs-membs (estimate-membs objs-feat-vals objs-groups initial-comms-props (range num-objs))
          _ (write-log "Initial object memberships found.")
          {est-comms-props :comms-props
           est-objs-membs :objs-membs} (run
                                         objs-feat-vals
                                         groups-objs
                                         objs-groups
                                         initial-objs-membs
                                         initial-comms-props
                                         max-iter)
          _ (write-log "Estimating all objects membership...")
          num-objs (mx/row-count objs-feat-vals)
          all-est-objs-membs (estimate-membs objs-feat-vals objs-groups est-comms-props (range num-objs))]
      
      (write-log "Saving output...")
      (save-matrix est-comms-props (format "%s/mle-comms-props.mm" out-dir))
      (save-matrix all-est-objs-membs (format "%s/mle-objs-membs.mm" out-dir)))))


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
