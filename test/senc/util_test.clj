(ns moc.util-test
  (:require [clojure.test :refer :all]
            [moc.util :refer :all]
            [clojure.core.matrix :as mx]
            [schema.test]))

(use-fixtures :once schema.test/validate-schemas)

(mx/set-current-implementation :vectorz)

(defn close? [delta-thresh expected actual]
  (< (Math/abs (- expected actual)) delta-thresh))

(defn close-seq? [delta-thresh expected actual]
  (every? (fn [[e a]] (< (Math/abs (- e a)) delta-thresh))
          (map vector expected actual)))
