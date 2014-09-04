(ns moc.experiment.optimize-maxlikelihood
  (:require [clojure.core.matrix :as mx]
            [schema.macros :as sm]
            [schema.core :as sc]
            [moc.objective :as mo])
  (:use [moc.schema :only [Vec Mat]])
  (:import (lbfgsb DifferentiableFunction
                   Minimizer
                   Bound
                   FunctionValues)))

;;; This was a poor attempt at using ML to estimate communities.
;;; Some of this code, such as the objective fns and gradient fns
;;; could be generalized and maybe useful for other purposes.

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

(sm/defn comm-objective-f :- DifferentiableFunction
  "Create objective function for single community."
  [comm-props :- Mat
   obj-mships :- Vec
   obj-feat :- Vec
   comm-idx :- Integer]
  (mo/obj-fn (partial mo/f obj-mships obj-feat comm-props comm-idx) 
             (conj (mo/f-dels-prop obj-mships comm-props comm-idx obj-feat)
                   mo/f-del-l-mult)))

(sm/defn comm-objective-h :- DifferentiableFunction
  "Create objective function for single community."
  [comm-props :- Mat
   obj-mships :- Vec
   obj-feat :- Vec
   comm-idx :- Integer]
  (mo/obj-fn (partial mo/h obj-mships obj-feat comm-props comm-idx) 
             (conj (mo/h-dels-prop obj-mships obj-feat)
                   mo/h-del-l-mult)))

(sm/defn comm-objective-h-approx-dels :- DifferentiableFunction
  "Create objective function for single community."
  [comm-props :- Mat
   obj-mships :- Vec
   obj-feat :- Vec
   comm-idx :- Number]
  (let [num-vars (inc (mx/dimension-count obj-feat 0))]
    (mo/obj-fn (partial mo/h obj-mships obj-feat comm-props comm-idx)
               (map #(partial mo/approx-del (partial mo/h obj-mships obj-feat comm-props comm-idx) %1 %2) (range num-vars) (repeat (double 0.00001))))))

(sm/defn bound :- Bound
  [[l u] :- [(sc/one (sc/maybe Double) "l") (sc/one (sc/maybe Double) "u")]]
  (Bound. l u))

(sm/defn comm-minimizer :- Minimizer
  [bound-pairs :- [[(sc/one (sc/maybe Double) "lower") (sc/one (sc/maybe Double) "upper")]]]
  (doto (Minimizer.)
    (.setBounds (map bound bound-pairs))))

(def of-f (comm-objective-f comm-props-mat (mx/get-row obj-memb-mat 0) (mx/get-row obj-feat-mat 0) (int 0)))
(def m-f (comm-minimizer (conj (vec (repeat num-dims [0.001 1.0])) [0.5 nil])))
;(.run m-f of-f (double-array (conj (into [] (mx/get-row comm-props-mat 0)) 0.001)))

(def of-h (comm-objective-h comm-props-mat (mx/get-row obj-memb-mat 0) (mx/get-row obj-feat-mat 0) (int 0)))
(def of-h (comm-objective-h truth-props-mat (mx/get-row obj-memb-mat 0) (mx/get-row obj-feat-mat 0) (int 0)))
(def m-h (comm-minimizer (conj (vec (repeat num-dims [0.001 1.0])) [0.001 nil])))

(def of-h-a (comm-objective-h-approx-dels comm-props-mat (mx/get-row obj-memb-mat 0) (mx/get-row obj-feat-mat 0) (int 0)))
(def m-h-a (comm-minimizer (conj (vec (repeat num-dims [(double 0.001) (double 1.0)])) 
                                 [(double 0.0001) nil])))

(defn try-it []
  (use '[moc.experiment.known-obj])
  (require '[clojure.core.matrix :as mx])
  (require '[schema.core :as sc])
  (sc/with-fn-validation 
    (.run m-h-a of-h-a (double-array (conj (into [] (mx/get-row comm-props-mat 0)) 10000)))))

; (sc/with-fn-validation (.run m-h of-h (double-array (conj (into [] (mx/get-row comm-props-mat 0)) 0.001))))

;; norm.vec(c(0, 0, 0, 0.075, 0.05, 0.05, 0.025, 0.025, 0, 0)),
;; norm.vec(c(0, 0, 0, 0.025, 0.050, 0.050, 0.075, 0.075, 0.100, 0.100)),
;; norm.vec(c(0.3, 0.2, 0.1, 0.05, 0.05, 0, 0, 0, 0, 0))),

;; Testing the moc obj
;(.getValues of (double-array (mx/get-row comm-props-mat 0)))

;; A test for l-bfgs-b wrapper
;; This works...
(comment (def q (proxy [DifferentiableFunction] []
                  (getValues [pt] 
                    (let [x (aget pt 0)]
                      (FunctionValues. (Math/pow (+ x 4) 2)
                                       (double-array [(* 2 (+ x 4))]))))))

         (def alg (doto (Minimizer.)
                    (.setBounds (map bound [[10.0 nil]]))))

         (def result (.run alg q (double-array [40]))))

;(def partials (conj (mo/f-dels-prop (mx/get-row obj-memb-mat 0) comm-props-mat 0 (mx/get-row obj-feat-mat 0)) mo/f-del-l-mult))
;(def grad-fn (apply juxt partials))
;(grad-fn (double-array (conj (into [] (mx/get-row comm-props-mat 0)) 0.001)))

;(def o-f-fn (partial mo/f (mx/get-row obj-memb-mat 0) (mx/get-row obj-feat-mat 0) comm-props-mat 0))
;(o-fn (double-array (conj (into [] (mx/get-row comm-props-mat 0)) 0.5)))

;(def o-h-fn (partial mo/h (mx/get-row obj-memb-mat 0) (mx/get-row obj-feat-mat 0) comm-props-mat 0))
;(o-h-fn (double-array (conj (into [] (mx/get-row comm-props-mat 0)) 0.5)))
