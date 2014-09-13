(ns moc.util
  (:require [clojure.core.matrix :as mx]
            [schema.core :as sc]
            [schema.macros :as sm])
  (:use [moc.schema :only [Vec Mat ProbVec BinVec]]))

(sm/defn proportional :- Vec
  "Normalize vector by L1-norm."
  [v :- Vec]
  (let [l1-norm (mx/ereduce (fn [sum x] (+ sum (Math/abs x))) v)]
    (if (zero? l1-norm)
      v
      (mx/div v l1-norm))))

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

  (loop [new-m (mx/sparse (apply mx/zero-matrix (mx/shape m)))
         row-idxs (mx/non-zero-indices row-indicator)]
    (if (empty? row-idxs)
      (mx/sparse new-m)
      (recur (mx/set-row new-m (first row-idxs) (mx/get-row m (first row-idxs)))
             (rest row-idxs)))))

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

