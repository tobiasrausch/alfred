bamStats installation
---------------------

`git clone --recursive https://github.com/tobiasrausch/bamStats.git`

`cd bamStats/`

`make all`

Running bamStats
------------------

`./src/bamStats -r <ref.fa> -o outprefix <align.bam>`

Plotting the insert size distribution.

`Rscript R/isize.R outprefix.isize.tsv`

Plotting the coverage distribution.

`Rscript R/coverage.R outprefix.coverage.tsv`

To convert all the alignment metrics from column format to rows to easily read it on screen.

`cat outprefix.metrics.tsv | datamash transpose | column -t`





