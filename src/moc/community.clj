(ns moc.community
  (:require [schema.core]
            [schema.macros :as sm])
  (:gen-class))

(sm/defrecord Obj
    [id :- Integer
     feat-counts :- [Double]])

(sm/defn obj
  [id :- Integer
   feat-counts :- [Double]]
  (Obj. id feat-counts))

(sm/defrecord Community
    [id :- Integer
     feat-props :- [Double]
     members :- [moc.community.Obj]])

(sm/defn community :- Community
  [id :-  Integer
   feat-props :- [Double]
   members :- [moc.community.Obj]]
  (Community. id feat-props members))
