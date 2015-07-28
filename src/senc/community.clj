(ns senc.community
  (:require [schema.core :as s]))

(s/defrecord Obj
    [id :- Integer
     feat-counts :- [Double]])

(s/defn obj
  [id :- Integer
   feat-counts :- [Double]]
  (Obj. id feat-counts))

(s/defrecord Community
    [id :- Integer
     feat-props :- [Double]
     members :- [moc.community.Obj]])

(s/defn community :- Community
  [id :-  Integer
   feat-props :- [Double]
   members :- [moc.community.Obj]]
  (Community. id feat-props members))
