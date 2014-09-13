(ns moc.estimate.mle
  (:require [clojure.core.matrix :as mx]
            [schema.core :as sc]
            [schema.macros :as sm])
  (:use [moc.schema :only [Vec Mat ProbVec]]
        [moc.util :only [proportional select-rows]]))

(sm/defn safe-log :- sc/Num
  "A log function that returns zero for input value zero.
  Simplifies working with probability vectors where zero-value
  elements can be ignored."
  [x :- sc/Num]
  (if (zero? x)
    0
    (Math/log x)))

(sm/defn solve-linear :- sc/Num
  "Solve a simple linear equation of the form y = ax + b."
  [y :- sc/Num
   a :- sc/Num
   b :- sc/Num]
  (/ (- y b) a))

(sm/defn estimate-mixed-mult :- ProbVec
  "We are solving a series of linear equations, 
  one for each categorical dimension.

  The observations, other community proportions, and 
  mixing weight of the community being estimated are 
  provided."
  [obs :- Vec
   other-props :- (sc/either ProbVec Vec)
   mix-weight :- sc/Num]
  (comment (mx/matrix (map (fn [y a b]
                             ;; If observations are already accounted for by
                             ;; other community distributions, then we know
                             ;; this community should include a near-zero
                             ;; probability for the dimension
                             ;; TODO: Is it not OK to just use 0?
                             (if (> y b)
                               (solve-linear y a b)
                               0))
                           (proportional obs) (repeat mix-weight) other-props)))

  ;; final proportional needed when there were negative values from subtracting the
  ;; proportion of observed feature values by the mixed other community parameters
  (->> (-> (proportional obs)
           (mx/sub other-props)
           (mx/div mix-weight))
       (map (partial max 0))
       (proportional)))

(sm/defn estimate-comm :- Vec
  [select-objs :- (sc/either (sm/=> [sc/Int] sc/Int) [[sc/Int]])
   obj-feats :- Mat
   obj-memb :- Mat
   comm-props :- Mat
   comm-idx :- sc/Int]
  (let [obj-idxs (select-objs comm-idx)
        comb-memb (proportional (reduce mx/add (select-rows obj-memb obj-idxs)))
        comb-feats (reduce mx/add (select-rows obj-feats obj-idxs))
        other-props (mx/mmul (mx/mset comb-memb comm-idx 0) comm-props)
        mix-weight (mx/mget comb-memb comm-idx)]
    (estimate-mixed-mult comb-feats other-props mix-weight)))

(sm/defn estimate-mix-from-terms :- Vec
  "Estimate the mixture of communities for an object."
  [obj-feats :- Mat
   obj-memb :- Mat
   comm-props :- Mat
   obj-idx :- sc/Int])
