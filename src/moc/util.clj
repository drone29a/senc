(ns moc.util
  (:require [clojure.core.matrix :as mx]
            [schema.core :as sc]
            [schema.macros :as sm])
  (:use [moc.schema :only [Vec Mat]]))

(sm/defn proportional :- Vec
  "Normalize vector by L1-norm."
  [v :- Vec]
  (let [l1-norm (mx/ereduce (fn [sum x] (+ sum (Math/abs x))) v)]
    (if (zero? l1-norm)
      v
      (mx/div v l1-norm))))

(sm/defn select-rows :- [Vec]
  [m :- Mat
   idxs :- [sc/Int]]
  (map (partial mx/get-row m) idxs))
