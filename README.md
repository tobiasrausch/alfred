<p align="center">
<img width="230" height="100" src="https://raw.githubusercontent.com/tobiasrausch/alfred/master/alfred.png">
</p>

[![Build Status](https://travis-ci.org/dellytools/delly.svg?branch=master)](https://travis-ci.org/tobiasrausch/alfred)


Alfred installation
---------------------

`git clone --recursive https://github.com/tobiasrausch/alfred.git`

`cd alfred/`

`make all`

Running Alfred
------------------

`./src/alfred -r <ref.fa> -o outprefix <align.bam>`

Plotting the insert size distribution.

`Rscript R/isize.R outprefix.isize.tsv`

Plotting the coverage distribution.

`Rscript R/coverage.R outprefix.coverage.tsv`

To convert all the alignment metrics from column format to rows to easily read it on screen.

`cat outprefix.metrics.tsv | datamash transpose | column -t`


Running Alfred using a bed file of target regions
---------------------------------------------------

If target regions are provided, Alfred computes the average coverage for each target and the on-target rate.

`./src/alfred -r <ref.fa> -b <targets.bed> -o outprefix <align.bam>`

Plotting the on-target rate.

`Rscript R/ontarget.R outprefix.ontarget.tsv`

Plotting the fraction of targets above a given coverage threshold.

`Rscript R/percov.R stats.bedcov.tsv`
