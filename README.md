bamStats installation
---------------------

`git clone --recursive https://github.com/tobiasrausch/bamStats.git`

`cd bamStats/`

`make all`

Running bamStats
------------------

`./src/bamStats -r <ref.fa> -o outprefix <align.bam>`

`Rscript R/isize.R outprefix.isize.tsv`

`Rscript R/coverage.R outprefix.coverage.tsv`

To read the metrics file on screen you may want to transpose rows and columns and layout columns:

`cat outprefix.metrics.tsv | datamash transpose | column -t`





