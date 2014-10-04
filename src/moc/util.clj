(ns moc.util
  (:require [clojure.core.matrix :as mx]
            [schema.core :as sc]
            [schema.macros :as sm])
  (:use [moc.schema :only [Vec Mat ProbVec BinVec]])
  (:import [mikera.matrixx.impl SparseRowMatrix]))

(set! *warn-on-reflection* true)

(comment (defprotocol PSparseRows
           "Protocol for selecting "
           (non-zero-rows [m] "gets rows with non-zero elements and returns a map of row index to row vector"))

         (extend-protocol PSparseRows
           SparseRowMatrix
           (non-zero-rows
             [m]
             (let [nzi (mx/non-zero-indices m)
                   num-inds (count nzi)]
               (loop [rows {}
                      i 0]
                 (if (<= num-inds i)
                   rows
                   (let [row-idx (nth nzi i)
                         row (.getRow m row-idx)]
                     (recur (assoc rows row-idx row)
                            (inc i)))))))))

(sm/defn proportional :- Vec
  "Normalize vector by L1-norm."
  [v :- Vec]
  (let [^ints nz-idxs (mx/non-zero-indices v)
        ;; TODO: non-zero-indices isn't guaranteed to return int[],
        ;;       but it does for the sparse data types and we want speeeeed
        num-idxs (alength nz-idxs)
        l1-norm (loop [idx (int 0)
                       sum (double 0.0)]
                  (if (> num-idxs idx)
                    (recur (inc idx)
                           (+ sum (Math/abs (double (mx/mget v (aget nz-idxs idx))))))
                    sum))]
    (if (zero? l1-norm)
      v
      (do (mx/div! v l1-norm)
          v))))

(sm/defn binary-vector :- BinVec
  "Create a binary vector with ones corresponding to given indices."
  [size :- sc/Int
   selected-idxs :- [sc/Int]]
  (let [v (mx/zero-vector size)]
    (doseq [idx selected-idxs]
      (mx/mset! v idx 1))
    v))

(sm/defn b-or :- BinVec
  "Binary OR operation for BinVecs."
  [& vecs :- [BinVec]]
  (let [max-size (apply max (map mx/row-count vecs))]
    (->> (mapcat mx/non-zero-indices vecs)
         (reduce conj #{})
         seq
         (binary-vector max-size))))

(sm/defn select-rows :- [Vec]
  [m :- Mat
   idxs :- [sc/Int]]
  (map (partial mx/get-row m) idxs))

(sm/defn selected-rows :- Mat
  "Return a matrix with unselected rows to zero."
  [m :- Mat
   row-indicator :- BinVec]

  ;; TODO: Curious when/if this is much/any slower
  (comment
    (mx/inner-product (mx/sparse (mx/diagonal-matrix row-indicator)) m))

  (let [[nrows ncols] (mx/shape m)
        new-m (SparseRowMatrix/create nrows ncols)]
    (doseq [row-idx (mx/non-zero-indices row-indicator)]
      (mx/set-row! new-m row-idx (mx/get-row m row-idx)))
    new-m))

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

(sm/defn safe-log :- sc/Num
  "A log function that returns a very small number for input value zero.
  Simplifies working with probability vectors where zero-value
  elements can be ignored."
  [x :- sc/Num]
  (if (zero? x)
    1e-10
    (Math/log x)))

(sm/defn log-multicat :- sc/Num
  "This is a bastardized not-quite-legit log PMF. It's essentially
  a multinomial without the multinomial coefficient. This is useful
  for computing relative change in likelihood for a distribution
  we're trying to estimate."
  [p :- Vec
   xs :- Vec]
  (mx/mget ^mikera.vectorz.AScalar (mx/inner-product xs (mx/emap safe-log p))))
