(ns moc.core-test
  (:require [clojure.test :refer :all]
            [moc.core :refer :all]
            [moc.util :refer [proportional]]
            [clojure.core.matrix :as mx]
            [schema.test]))

(use-fixtures :once schema.test/validate-schemas)

(mx/set-current-implementation :vectorz)

(defn close? [delta-thresh expected actual]
  (< (Math/abs (- expected actual)) delta-thresh))

(defn close-seq? [delta-thresh expected actual]
  (every? (fn [[e a]] (< (Math/abs (- e a)) delta-thresh))
          (map vector expected actual)))

(deftest estimate-memb-test
  (is (close-seq? 1e-4
                  [0.3333 0.6666]
                  (estimate-memb (mx/sparse [10 0 20 0 30])
                                 (mx/sparse [[0.6 0.3 0.1 0 0]
                                             [0 0 0.1 0.3 0.6]])))))

(deftest restricted-comms-test
  (is (mx/equals (mx/matrix [[0.5 0.3 0.2]
                             [0.7 0.1 0.2]
                             [0 0 0]])
                 (restricted-comms (mx/matrix [[1 0 0]
                                               [1 1 0]
                                               [1 0 0]])
                                   (mx/matrix [[0.5 0.3 0.2]
                                               [0.7 0.1 0.2]
                                               [0.4 0.4 0.2]])
                                   #{0 1}))))

(deftest e-step-test
  (testing "trivial case of single membership"
    (is (mx/equals (mx/matrix [[1 0 0]
                               [0 1 0]
                               [0 0 1]])
                   (e-step (partial restricted-comms (mx/matrix [[1 0 0]
                                                                 [0 1 0]
                                                                 [0 0 1]]))
                           (mx/matrix [[5 2 0 0]
                                       [1 3 0 0]
                                       [0 0 3 10]])
                           (mx/matrix [[0.4 0.2 0.2 0.2]
                                       [0.2 0.4 0.2 0.0]
                                       [0.0 0.0 0.3 0.7]])))))
  (testing "mixture membership"
    (is (close-seq? 1e-2
                    (concat [0.612 0.388 0.0]
                            [0.433 0.566 0.0]
                            [0.521 0.0 0.479]
                            [0.0 0.41 0.59])
                    (mx/eseq (e-step (partial restricted-comms (mx/matrix [[1 1 0]
                                                                           [1 1 0]
                                                                           [1 0 1]
                                                                           [0 1 1]]))
                                     (mx/matrix [[5 2 1 1]
                                                 [1 3 1 0]
                                                 [5 2 3 10]
                                                 [5 2 3 10]])
                                     (mx/matrix [[0.4 0.2 0.2 0.2]
                                                 [0.2 0.4 0.2 0.0]
                                                 [0.0 0.0 0.3 0.7]])))))))

(deftest subgraph-membership-test
  (testing "trivial case"
    (is (close-seq? 1e-2
                    [0.84 0.12 0.04]
                    (subgraph-membership (mx/matrix [[1.0 0 0]
                                                     [0.6 0.3 0.1]])
                                         (mx/matrix [0.6
                                                     0.4]))))))

(deftest other-mixed-props-test
  (testing "trivial case when objects only belong to community being estimated"
    ;; For comm-idx 0
    (is (close-seq? 1e-2
                    [0.0 0.0 0.0]
                    (other-mixed-props (mx/matrix [0.0 0.0 0.0])
                                       (mx/matrix [[0.6 0.4 0.0]
                                                   [0.3 0.7 0.0]
                                                   [0.9 0.0 0.1]])))))
  (testing "simple case"
    ;; For comm-idx 0
    (is (close-seq? 1e-2
                    [0.30 0.28 0.02]
                    (other-mixed-props (mx/matrix [0.0 0.4 0.2])
                                       (mx/matrix [[0.6 0.4 0.0]
                                                   [0.3 0.7 0.0]
                                                   [0.9 0.0 0.1]]))))
    ;; For comm-idx 1
    (is (close-seq? 1e-2
                    [0.42 0.16 0.02]
                    (other-mixed-props (mx/matrix [0.4 0.0 0.2])
                                       (mx/matrix [[0.6 0.4 0.0]
                                                   [0.3 0.7 0.0]
                                                   [0.9 0.0 0.1]]))))))

(deftest estimate-comm-props-test
  (testing "trivial case of estimating community parameters"
    (is (close-seq? 1e-2
                    [0.66 0.33 0 0]
                    (mx/eseq (estimate-comm-props (mx/matrix [[1 0]
                                                              [0 1]])
                                                  (mx/matrix [[1.0 0]
                                                              [0 1.0]])
                                                  (mx/matrix [[10 5 0 0]
                                                              [0 0 10 10]])
                                                  (mx/matrix [[0.7 0.3 0 0]
                                                              [0 0 0.6 0.4]])
                                                  0)))))
  (testing "case of estimating community parameters where result sums greater 
            than one b/c some features occurred less than expected and 
            not were overaccoutned for by other community distributions"
    (is (close-seq? 1e-2
                    [0.69 0.31 0.0 0.0]
                    (mx/eseq (estimate-comm-props (mx/matrix [[1 0 0]
                                                              [0 0 1]
                                                              [1 0 1]])
                                                  (mx/matrix [[1.0 0 0]
                                                              [0 1.0 0]
                                                              [0.3 0 0.7]])
                                                  (mx/matrix [[10 5 0 0]
                                                              [0 0 10 10]
                                                              [6 2 4 10]])
                                                  (mx/matrix [[0.7 0.3 0 0]
                                                              [0 0 0.6 0.4]
                                                              [0 0 0.3 0.7]])
                                                  0))))))

;; num object groups = 5
;; num comms = 5
;; num features = 6
(deftest ^:slow run-uniform-test
  (testing "simple set of communities and objects"
    (let [feat-vals (mx/matrix [[73  8 19  0  0]
                                [72  1 27  0  0]
                                [69  1 30  0  0]
                                [65  9 26  0  0]
                                [68 22 10  0  0]
                                [72 14 14  0  0]
                                [68  0 32  0  0]
                                [74 16 10  0  0]
                                [66  5 29  0  0]
                                [64 22 14  0  0]
                                [71 11 18  0  0]
                                [75  5 20  0  0]
                                [66 17 17  0  0]
                                [76 13 11  0  0]
                                [66  1 33  0  0]
                                [77  7 16  0  0]
                                [70 11 19  0  0]
                                [72 14 14  0  0]
                                [69 14 17  0  0]
                                [66 12 22  0  0]
                                [68 13 19  0  0]
                                [68  0 32  0  0]
                                [67 10 23  0  0]
                                [70  6 24  0  0]
                                [76  1 23  0  0]
                                [67 19 14  0  0]
                                [64 16 20  0  0]
                                [62 22 16  0  0]
                                [74 12 14  0  0]
                                [65  0 35  0  0]
                                [65  5 30  0  0]
                                [65 15 20  0  0]
                                [71 16 13  0  0]
                                [68 15 17  0  0]
                                [76  1 23  0  0]
                                [68  9 23  0  0]
                                [65 21 14  0  0]
                                [70  6 24  0  0]
                                [72  8 20  0  0]
                                [70  0 30  0  0]
                                [80  6 14  0  0]
                                [77 13 10  0  0]
                                [74  6 20  0  0]
                                [73  9 18  0  0]
                                [63  5 32  0  0]
                                [76  0 24  0  0]
                                [67 16 17  0  0]
                                [66 21 13  0  0]
                                [76  0 24  0  0]
                                [67 14 19  0  0]
                                [71  8 21  0  0]
                                [71 18 11  0  0]
                                [71 11 18  0  0]
                                [70 22  8  0  0]
                                [78 11 11  0  0]
                                [75  6 19  0  0]
                                [71 13 16  0  0]
                                [72 16 12  0  0]
                                [73 15 12  0  0]
                                [75  6 19  0  0]
                                [32  3 20 17 28]
                                [66  4 25  2  3]
                                [56 11 17  6 10]
                                [56 11 15  8 10]
                                [43  3 21 11 22]
                                [27  4 17 23 29]
                                [55  4 21  9 11]
                                [50  4 24  6 16]
                                [55 17 25  0  3]
                                [45  5 20 13 17]
                                [75  9 16  0  0]
                                [51  7 17 11 14]
                                [37  2 31 14 16]
                                [64 10 18  3  5]
                                [64  0 22  3 11]
                                [55  1 27  7 10]
                                [52  6 15 14 13]
                                [61  7 21  5  6]
                                [53  5 22  4 16]
                                [37  4 20 13 26]
                                [63 14 11  8  4]
                                [75 13 11  1  0]
                                [63  1 29  3  4]
                                [69 14 16  0  1]
                                [34  0 24 15 27]
                                [47  8 17 11 17]
                                [63 14 12  4  7]
                                [68  5 26  1  0]
                                [49 11 16  8 16]
                                [68  4 22  2  4]
                                [37  0 21 10 32]
                                [61  0 28  7  4]
                                [60  0 30  2  8]
                                [57  0 30  3 10]
                                [47  0 28 13 12]
                                [70  0 29  0  1]
                                [63  0 36  0  1]
                                [49  0 21 12 18]
                                [51  0 26  4 19]
                                [71  0 29  0  0]
                                [51 32 14  2  1]
                                [48 24 18  7  3]
                                [45 24 16  7  8]
                                [34 57  5  0  4]
                                [43 25 11 15  6]
                                [61  8 20  8  3]
                                [27 10 21 28 14]
                                [54 16 25  3  2]])
          obj-groups (mx/matrix (concat (repeat 60 [1 1 0 0 0])
                                        (repeat 30 [1 1 1 0 0])
                                        (repeat 10 [0 1 1 0 0])
                                        (repeat 5 [1 0 0 1 1])
                                        (repeat 3 [0 1 0 1 1])))
          cores (mx/matrix [(concat (repeat 60 1) (repeat 48 0))
                            (concat (repeat 60 0) (repeat 30 1) (repeat 18 0))
                            (concat (repeat 90 0) (repeat 10 1) (repeat 8 0))
                            (concat (repeat 100 0) (repeat 5 1) (repeat 3 0))
                            (concat (repeat 105 0) (repeat 3 1))])
          init-comm-props (estimate-props feat-vals cores)
          init-obj-membs (mx/matrix (map (fn [obj-feat-vals]
                                           (estimate-memb obj-feat-vals init-comm-props))
                                         feat-vals))]
        (is (close-seq? 1e-3
             [0 0 0 0 0]
             (-> (run
                   feat-vals 
                   obj-groups
                   init-obj-membs
                   init-comm-props
                   50)
                 :comm-props
                 mx/eseq))))))

(deftest ^:slow run-select-heavy-test
  (testing "simple set of communities and objects"
    (let [feat-vals (mx/matrix [[76 15 9  0  0]   
                                [66 23 11 0  0]   
                                [70 21 9  0  0]   
                                [77 12 11 0  0]   
                                [73 13 14 0  0]   
                                [70 15 15 0  0]   
                                [77 11 12 0  0]   
                                [62 19 19 0  0]   
                                [71 19 10 0  0]   
                                [70 21 9  0  0]   
                                [66 20 14 0  0]   
                                [74 17 9  0  0]   
                                [70 23 7  0  0]   
                                [79 10 11 0  0]   
                                [71 14 15 0  0]   
                                [71 19 10 0  0]   
                                [70 20 10 0  0]   
                                [66 24 10 0  0]   
                                [81 9  10 0  0]   
                                [70 17 13 0  0]   
                                [70 18 12 0  0]   
                                [68 19 13 0  0]   
                                [71 17 12 0  0]   
                                [61 23 16 0  0]   
                                [68 24 8  0  0]   
                                [73 8  19 0  0]   
                                [70 19 11 0  0]   
                                [72 17 11 0  0]   
                                [68 25 7  0  0]   
                                [75 21 4  0  0]   
                                [68 23 9  0  0]   
                                [74 17 9  0  0]   
                                [70 21 9  0  0]   
                                [67 15 18 0  0]   
                                [73 18 9  0  0]   
                                [71 20 9  0  0]   
                                [71 17 12 0  0]   
                                [71 20 9  0  0]   
                                [70 17 13 0  0]   
                                [74 8  18 0  0]   
                                [74 14 12 0  0]   
                                [66 24 10 0  0]   
                                [71 19 10 0  0]   
                                [68 23 9  0  0]   
                                [62 23 15 0  0]   
                                [69 18 13 0  0]   
                                [71 16 13 0  0]   
                                [74 14 12 0  0]   
                                [79 13 8  0  0]   
                                [70 25 5  0  0]   
                                [71 19 10 0  0]   
                                [75 14 11 0  0]   
                                [75 15 10 0  0]   
                                [72 16 12 0  0]   
                                [61 24 15 0  0]   
                                [67 28 5  0  0]   
                                [59 23 18 0  0]   
                                [67 19 14 0  0]   
                                [73 16 11 0  0]   
                                [75 14 11 0  0]   
                                [60 4  24 3  9]   
                                [65 1  32 1  1]   
                                [66 3  30 0  1]   
                                [71 0  28 1  0]   
                                [60 1  38 0  1]   
                                [72 0  21 4  3]   
                                [76 0  21 1  2]   
                                [63 1  31 4  1]   
                                [70 0  26 3  1]   
                                [67 0  26 5  2]   
                                [70 1  26 1  2]   
                                [74 0  25 0  1]   
                                [54 4  28 5  9]   
                                [60 2  33 3  2]   
                                [70 1  28 0  1]   
                                [68 2  27 1  2]   
                                [66 5  26 1  2]   
                                [77 1  21 1  0]   
                                [63 1  27 4  5]   
                                [66 0  31 1  2]   
                                [71 0  21 5  3]   
                                [61 0  35 1  3]   
                                [64 1  30 2  3]   
                                [69 1  29 1  0]   
                                [74 8  17 0  1]   
                                [64 1  33 1  1]   
                                [63 4  29 2  2]   
                                [69 0  30 1  0]   
                                [63 1  36 0  0]   
                                [70 6  23 1  0]   
                                [33 0  21 15 31]  
                                [35 0  21 21 23]  
                                [21 0  27 21 31]  
                                [40 0  22 7  31]  
                                [33 0  30 15 22]  
                                [27 0  25 20 28]  
                                [28 0  29 13 30]  
                                [28 0  26 18 28]  
                                [39 0  21 10 30]  
                                [39 0  18 18 25]  
                                [36 42 16 4  2]   
                                [40 54 5  0  1]   
                                [35 44 15 3  3]   
                                [32 56 9  2  1]   
                                [35 40 16 7  2]   
                                [22 16 25 27 10]  
                                [23 15 15 35 12]  
                                [13 19 21 26 21]])
          obj-groups (mx/matrix (concat (repeat 60 [1 1 0 0 0])
                                        (repeat 30 [1 1 1 0 0])
                                        (repeat 10 [0 1 1 0 0])
                                        (repeat 5 [1 0 0 1 1])
                                        (repeat 3 [0 1 0 1 1])))
          cores (mx/matrix [(concat (repeat 60 1) (repeat 48 0))
                            (concat (repeat 60 0) (repeat 30 1) (repeat 18 0))
                            (concat (repeat 90 0) (repeat 10 1) (repeat 8 0))
                            (concat (repeat 100 0) (repeat 5 1) (repeat 3 0))
                            (concat (repeat 105 0) (repeat 3 1))])
          init-comm-props (estimate-props feat-vals cores)
          init-obj-membs (mx/matrix (map (fn [obj-feat-vals]
                                           (estimate-memb obj-feat-vals init-comm-props))
                                         feat-vals))]
        (is (close-seq? 1e-3
             [0 0 0 0 0]
             (-> (run
                   feat-vals 
                   obj-groups
                   init-obj-membs
                   init-comm-props
                   50)
                 :comm-props
                 mx/eseq))))))
