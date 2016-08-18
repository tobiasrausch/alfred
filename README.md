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
	


