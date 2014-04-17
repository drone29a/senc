(ns moc.schema
  (:require [clojure.core.matrix :as mx]
            [schema.core :as sc]))

(def Vec
  "Any vector"
  (sc/pred mx/vec? 'vec?))

(def Mat
  "Any matrix"
  (sc/pred mx/matrix? 'matrix?))
