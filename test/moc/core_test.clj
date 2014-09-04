(use-fixtures :once schema.test/validate-schemas)

(ns moc.core-test
  (:require [clojure.test :refer :all]
            [moc.core :refer :all]
            [matrix.core :as mx]))

(defn close? [delta-thresh expected actual]
  (< (Math/abs (- expected actual)) delta-thresh))

(defn close-seq? [delta-thresh expected actual]
  (every? (fn [[e a]] (< (Math/abs (- e a)) delta-thresh))
          (map vector expected actual)))

(deftest normalize-log-prob-test
  (is (close-seq? 1e-3
                  [0.4223,0.4223,0.1553]
                  (mx/eseq (normalize-log-prob (mx/matrix [-111 -111 -112]))))))
