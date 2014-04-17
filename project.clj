(defproject moc "0.1.0-SNAPSHOT"
  :description "Modelling communities and topics in networks."
  :url "http://example.com/FIXME"
  :license {:name "MIT"
            :url "http://opensource.org/licenses/MIT"}
  :dependencies [[org.clojure/clojure "1.5.1"]
                 [org.clojure/tools.cli "0.2.1"]
                 [prismatic/schema "0.2.1"]
                 [notion "0.1.0-SNAPSHOT"]
                 [munge "0.1.0-SNAPSHOT"]
                 [net.mikera/core.matrix "0.20.1-SNAPSHOT"]
                 [net.mikera/vectorz-clj "0.20.1-SNAPSHOT"]
                 [aysylu/loom "0.4.2"]
                 [mkobos/lbfgsb-wrapper "1.1.3-SNAPSHOT" 
                  :native-prefix ""]]
  :native-path "native"
  :jvm-opts ["-Xmx4g" "-server" "-XX:+UseConcMarkSweepGC" "-Djava.library.path=native/"]
  :global-vars {*warn-on-reflection* true})
