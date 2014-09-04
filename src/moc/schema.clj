(ns moc.schema
  (:require [clojure.core.matrix :as mx]
            [schema.core :as sc]))

(def Vec
  "Any vector"
  (sc/pred mx/vec? 'vec?))

(def Mat
  "Any matrix"
  (sc/pred mx/matrix? 'matrix?))

(def ProbVec
  "Probabilty vector"
  (sc/pred (fn [x] (and (mx/vec? x)
                        (< 0.001
                           (Math/abs (- 1 (mx/esum x))))))
           'prob-vec?))
