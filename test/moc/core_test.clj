(use-fixtures :once schema.test/validate-schemas)

(ns moc.core-test
  (:require [clojure.test :refer :all]
            [moc.core :refer :all]))

(deftest a-test
  (testing "FIXME, I fail."
    (is (= 0 1))))
