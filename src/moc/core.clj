(ns moc.core
  (:use [clojure.tools.cli :only [cli]])
  (:gen-class))

(defn load-data [adj-path obj-feat-path cores-path]
  )

(defn run
  "Run the algorithm. Main calls this, so can you."
  [adj-mat obj-feat-mat cores]
  )

(defn -main
  [& args]
  (let [[opts args banner] (cli args
                                ["-h" "--help" "Print this help message and exit" :default false :flag true]
                                ["-o" "--out-dir" "Output directory" :default (System/getProperty "user.dir")])]
    (when (or (:help opts) (not= 3 (count args)))
      (println banner)
      (System/exit 0))

    (let [adj-path (nth args 0) ; object links in CTF
          obj-feat-path (nth args 1) ; object x (real-valued) feature
          cores-path (nth args 2) ; csv of obj IDs, one core per line
          [adj-mat obj-feat-mat cores] (load-data adj-path obj-feat-path cores-path)]    
      
      ;; Return cluster ids
      {:clusts (run adj-mat obj-feat-mat cores)})))
