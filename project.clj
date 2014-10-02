(defproject moc "0.1.0-SNAPSHOT"
  :description "Modelling communities and topics in networks."
  :url "http://example.com/FIXME"
  :license {:name "MIT"
            :url "http://opensource.org/licenses/MIT"}
  :dependencies [[org.clojure/clojure "1.6.0"]
                 [org.clojure/tools.cli "0.3.1"]
                 [prismatic/schema "0.2.6"]
                 [notion "0.1.0-SNAPSHOT"]
                 [munge "0.1.0-SNAPSHOT"]
                 [net.mikera/core.matrix "0.30.2-SNAPSHOT"]
                 [net.mikera/vectorz-clj "0.25.1-SNAPSHOT"]
                 [aysylu/loom "0.5.1-SNAPSHOT"]
                 [mkobos/lbfgsb-wrapper "1.1.3-SNAPSHOT" 
                  :native-prefix ""]]
  :native-path "native"
  :main moc.core
  :jvm-opts ["-agentpath:/Applications/YourKit_Java_Profiler_2014_build_14104.app/bin/mac/libyjpagent.jnilib"
             "-Xdebug" "-Xmx8g" "-server" "-XX:+UseConcMarkSweepGC" "-Djava.library.path=native/"]
;;  :jvm-opts ["-Xmx8g" "-server" "-XX:+UseConcMarkSweepGC" "-Djava.library.path=native/"]
  :global-vars {*warn-on-reflection* true}
  :test-selectors {:default (complement :slow)
                   :slow :slow
                   :all (constantly true)})
