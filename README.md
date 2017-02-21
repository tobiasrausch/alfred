alfred installation
---------------------

`git clone --recursive https://github.com/tobiasrausch/alfred.git`

`cd alfred/`

`make all`

Running alfred
------------------

`./src/alfred -r <ref.fa> -o outprefix <align.bam>`

Plotting the insert size distribution.

`Rscript R/isize.R outprefix.isize.tsv`

Plotting the coverage distribution.

`Rscript R/coverage.R outprefix.coverage.tsv`

To convert all the alignment metrics from column format to rows to easily read it on screen.

`cat outprefix.metrics.tsv | datamash transpose | column -t`


Running alfred using a bed file of target regions
---------------------------------------------------

If target regions are provided, alfred computes the average coverage for each target and the on-target rate.

`./src/alfred -r <ref.fa> -b <targets.bed> -o outprefix <align.bam>`

Plotting the on-target rate.

`Rscript R/ontarget.R outprefix.ontarget.tsv`

Plotting the fraction of targets above a given coverage threshold.

`Rscript R/percov.R stats.bedcov.tsv`
