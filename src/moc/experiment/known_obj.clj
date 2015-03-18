(ns moc.experiment.known-obj
  (:require [clojure.core.matrix :as mx]
            [schema.macros :as sm]
            [schema.core :as sc]
            [moc.objective :as mo])
  (:use [moc.schema :only [Vec Mat]])
  (:import (lbfgsb DifferentiableFunction
                   Minimizer
                   Bound
                   FunctionValues)))

;;;; The object membership proportion vectors are known
;;;; and we attempt to estimate the categorical distribution
;;;; associated with each community.

(mx/set-current-implementation :vectorz)

(sm/defn proportional :- Vec
  "Normalize vector by L1-norm."
  [v :- Vec]
  (let [l1-norm (mx/ereduce (fn [sum x] (+ sum (Math/abs x))) v)]
    (if (zero? l1-norm)
      v
      (mx/div v l1-norm))))

(def num-dims 10) ; k - num of categorical dimensions
(def num-objs 3) ; n - num of objects
(def num-comms 3) ; m - num of communities

(def obj-adj-mat (mx/matrix [[0 1 1]
                             [1 0 1]
                             [1 1 0]]))
(def obj-feat-mat (mx/matrix [[0 0 0 147 131 142 99 106 52 58]
                              [77 37 24 89 75 62 30 30 0 0]
                              [123 91 44 35 68 24 40 32 45 56]]))
(def obj-memb-mat (mx/matrix [[0.5 0.5 0.0001]
                              [0.5 0.0001 0.5]
                              [0.0001 0.5 0.5]]))

(def truth-props-mat (mx/matrix [(proportional [0 0 0 0.075 0.05 0.05 0.025 0.025 0 0])
                                 (proportional [0 0 0 0.025 0.050 0.050 0.075 0.075 0.100 0.100])
                                 (proportional [0.3 0.2 0.1 0.05 0.05 0 0 0 0 0])]))

(sm/defn weighted-prop :- Vec
  [weights :- Vec
   counts :- Mat]
  (let [weighted-counts (mx/mmul weights counts)]
    (mx/div weighted-counts (mx/esum weighted-counts))))

(sm/defn estimate-comms :- Mat
  [obj-feat-mat :- Mat
   obj-memb-mat :- Mat]
  (mx/matrix (map weighted-prop 
                  (mx/columns obj-memb-mat) 
                  (repeat obj-feat-mat))))

(def comm-props-mat (estimate-comms obj-feat-mat obj-memb-mat))

(comment
  (defn try-it []
    (use '[moc.experiment.known-obj])
    (require '[clojure.core.matrix :as mx])
    (require '[schema.core :as sc])
    (sc/with-fn-validation 
      (.run m-h-a of-h-a (double-array (conj (into [] (mx/get-row comm-props-mat 0)) 10000))))))

;;; scratch

(comment

  (def data (load-data "/Users/matt/data/cycle_comms/1_objs-affinity.csv"
                       "/Users/matt/data/cycle_comms/1_objs-feats.mtx"
                       "/Users/matt/data/cycle_comms/1_objs-membership.mtx"
                       "/Users/matt/data/cycle_comms/1_comm-params.mtx")))
