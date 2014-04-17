(ns moc.core
  (:require [clojure.core.matrix :as mx]
            [schema.core :as sc]
            [schema.macros :as sm])
  (:use [clojure.tools.cli :only [cli]]
        [moc.community :only [community obj]]
        [munge.io.matrix-ctf :only [load-matrix]]
        [moc.schema :only [Vec Mat]])
  (:gen-class))

#_(set-current-implementation :vectorz)

(defn load-data 
  [obj-adj-path 
   obj-feat-path 
   cores-path]
  [(load-matrix obj-adj-path)
   (load-matrix obj-feat-path)
   (load-matrix cores-path)])

;;; Need a way to estimate obj membership weights to each community/core
;obj-memb-mat <- (obj-memb cores obj-adj-mat)

(sm/defn partial-delta :- Number 
  "Calculate the partial derivative of the log-likelihood
  for the given categorical dimension."
  [memb-weight :- Number
   other-memb-weights :- Vec
   comm-prop :- Number
   all-comms-prop :- Vec
   feat-count :- Number]
  (let [denom (+ (* memb-weight comm-prop)
                 (mx/dot other-memb-weights all-comms-prop))
        denom (if (zero? denom) 0.0001 denom)]
    (* feat-count
       (/ memb-weight denom))))

(sm/defn gradient :- Vec
  "Caclulates the gradient for a community distribution. This is done as likelihood function
  of the mixed multinomial for a single object."
  [obj-feat-vec :- Vec
   obj-memb-vec :- Vec
   comm-props-mat :- Mat
   comm-idx :- sc/Int]
  (let [memb-weight (mx/mget obj-memb-vec comm-idx)
        other-memb-weights (mx/mset obj-memb-vec comm-idx 0)
        ptl-delta (partial partial-delta memb-weight other-memb-weights)]
    (vec (map (fn [dim-idx]
                (let [feat-count (mx/mget obj-feat-vec dim-idx)
                      comm-prop (mx/mget comm-props-mat comm-idx dim-idx)
                      all-comms-prop (mx/get-column comm-props-mat dim-idx)]
                  (ptl-delta comm-prop
                             all-comms-prop
                             feat-count)))
              (range (mx/column-count comm-props-mat))))))

(sm/defn proportional :- Vec
  "Normalize vector by L1-norm."
  [v :- Vec]
  (let [l1-norm (mx/ereduce (fn [sum x] (+ sum (Math/abs x))) v)]
    (if (zero? l1-norm)
      v
      (mx/div v l1-norm))))

(sm/defn proportional! :- Vec
  "Normalize vector by L1-norm."
  [v :- Vec]
  (let [l1-norm (mx/ereduce (fn [sum x] (+ sum (Math/abs x))) v)]
    (if (zero? l1-norm)
      v
      (mx/div! v l1-norm))))

(sm/defn select-obj :- sc/Int
  "Select an object for optimizing a community.
  For now, we only consider the membership weights of all objects
  for a single community."
  [objs-vec :- Vec]
  #_(rand (mx/row-count objs-vec))
  (first (reduce (fn [[max-idx max-val] [next-idx next-val]]
                   (if (> next-val max-val)
                     [next-idx next-val]
                     [max-idx max-val]))
                 ;; The index-seq returns indices in a vec, so we conj.
                 ;; Shuffle so we have a chance of gettin different obj when multiple
                 ;; have same mixing.
                 (shuffle (map conj (mx/index-seq objs-vec) (mx/eseq objs-vec))))))

(sm/defn run
  "Run the algorithm. Main calls this, so can you."
  [obj-adj-mat  :- Mat
   obj-feat-mat :- Mat
   obj-memb-mat :- Mat
   comm-props-mat :- Mat
   max-iter :- sc/Int
   update-rate :- Number]

  (loop [memb obj-memb-mat
         props comm-props-mat
         iter-count 0]
    (if (= max-iter iter-count)      
      comm-props-mat
      (let [obj-idxs (map select-obj (mx/columns obj-memb-mat))
            num-comms (mx/row-count comm-props-mat)
            updates (doall (map gradient 
                                (map (partial mx/get-row obj-feat-mat) obj-idxs)
                                (map (partial mx/get-row obj-memb-mat) obj-idxs)
                                (repeat comm-props-mat)
                                (range num-comms)))]
        ;;; Need to be sure all updates have been calculated before updating
        ;;; since we're currently mutating. Maybe make a copy to work with? 
        ;;; Could be immutable copy.
        (doseq [i (-> comm-props-mat mx/row-count range)]
          ;; (mx/set-row! props i (mx/+ (mx/get-row props i)
          ;;                            (nth updates i)))
          #_(println (mx/scale (nth updates i) update-rate))
          (doto (mx/get-row props i)
            ;; NOTE: We are adding the updates since we are maximizing/ascending
            (mx/add! (mx/emul (nth updates i) update-rate))
            proportional!))
        (recur memb props (inc iter-count))))))

(defn -main
  [& args]
  (let [[opts args banner] (cli args
                                ["-h" "--help" "Print this help message and exit" :default false :flag true]
                                ["-o" "--out-dir" "Output directory" :default (System/getProperty "user.dir")])]
    (when (or (:help opts) (not= 3 (count args)))
      (println banner)
      (System/exit 0))

    (let [obj-adj-path (nth args 0) ; object links in CTF
          obj-feat-path (nth args 1) ; object x (real-valued) feature
          cores-path (nth args 2) ; csv of obj IDs, one core per line
          [obj-adj-mat obj-feat-mat cores] (load-data obj-adj-path obj-feat-path cores-path)]    
      
      ;; Return cluster ids
      {:clusts (run obj-adj-mat obj-feat-mat cores)})))
