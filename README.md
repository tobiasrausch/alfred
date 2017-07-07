<p align="center">
<img width="230" height="100" src="https://raw.githubusercontent.com/tobiasrausch/alfred/master/alfred.png">
</p>

[![Build Status](https://travis-ci.org/dellytools/delly.svg?branch=master)](https://travis-ci.org/tobiasrausch/alfred)
[![Docker Automated buil](https://img.shields.io/docker/automated/jrottenberg/ffmpeg.svg?style=flat-square)](https://hub.docker.com/r/trausch/alfred/)


Alfred installation
---------------------

`git clone --recursive https://github.com/tobiasrausch/alfred.git`

`cd alfred/`

`make all`

Running Alfred
------------------

`./src/alfred -r <ref.fa> -o outprefix <align.bam>`

Plotting basic quality control metrics
--------------------------------------

Plotting the mean base quality across the read

`Rscript R/basequal.R outprefix.basequal.tsv`

Plotting the base composition across the read

`Rscript R/contentACGTN.R outprefix.contentACGTN.tsv`

Plotting the read length distribution

`Rscript R/readlength.R outprefix.readlength.tsv`

Plotting the mapping quality distribution

`Rscript R/mapq.R outprefix.mapq.tsv`

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

For instance, for a human whole-exome data set.

`./src/alfred -r <hg19.fa> -b exon/exon.hg19.bed.gz -o outprefix <exome.bam>`

Plotting the on-target rate.

`Rscript R/ontarget.R outprefix.ontarget.tsv`

Plotting the fraction of targets above a given coverage threshold.

`Rscript R/bedcov.R outprefix.bedcov.tsv`


Alfred wrapper scripts
----------------------

Alfred contains a wrapper script to combine all QC plots and metrics into a single pdf.

`./alfred.sh <ref.fa> <outprefix> <align.bam>`

A bed file of target regions is optional.

`./alfred.sh <ref.fa> <outprefix> <align.bam> <exome.bed>`


Example E. coli data set
------------------------

`./alfred.sh exampledata/E.coli.fa.gz exampledata/out exampledata/E.coli.cram`

The final pdf is exampledata/out.pdf


Example plots
-------------


[Whole-genome paired-end data with multiple read groups](https://raw.githubusercontent.com/tobiasrausch/alfred/master/exampleplots/NA06985.pe.pdf)

[Whole-exome paired-end data with on-target rate](https://raw.githubusercontent.com/tobiasrausch/alfred/master/exampleplots/HG00112.wes.pdf)

[Whole-genome jumping library](https://raw.githubusercontent.com/tobiasrausch/alfred/master/exampleplots/HG00513.mp.pdf)
