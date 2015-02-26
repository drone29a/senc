(defproject moc "0.1.0-SNAPSHOT"
  :description "Modelling communities and topics in networks."
  :url "http://example.com/FIXME"
  :license {:name "MIT"
            :url "http://opensource.org/licenses/MIT"}
  :dependencies [[org.clojure/clojure "1.7.0-alpha5"]
                 [org.clojure/tools.cli "0.3.1"]
                 [prismatic/schema "0.3.3"]
                 [notion "0.1.0-SNAPSHOT"]
                 [munge "0.1.0-SNAPSHOT"]
                 ;; [net.mikera/core.matrix "0.32.2-SNAPSHOT"]
                 ;; [net.mikera/vectorz-clj "0.28.1-SNAPSHOT"]
                 [net.mikera/core.matrix "0.33.3-SNAPSHOT"]
                 [net.mikera/vectorz-clj "0.29.1-SNAPSHOT"]
                 [com.climate/claypoole "0.4.0"]
                 ;; TODO: remove once objective ns moved to notion.
                 [mkobos/lbfgsb-wrapper "1.1.3-SNAPSHOT" :native-prefix ""]
                 ]
  :native-path "native"
  :main moc.core
  :profiles {:common {:jvm-opts ["-Xmx14g" "-server" "-XX:+UseConcMarkSweepGC" "-Djava.library.path=native/" "-XX:+TieredCompilation"]}
             :dev [:common {:jvm-opts ["-agentpath:/Applications/YourKit_Java_Profiler_2014_build_14116.app/Contents/Resources/bin/mac/libyjpagent.jnilib" "-Xdebug"]}]
             :prod [:common]}
  :global-vars {*warn-on-reflection* true}
  :test-selectors {:default (complement :slow)
                   :slow :slow
                   :all (constantly true)})
