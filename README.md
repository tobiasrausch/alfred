<p align="center">
   <img width="450" src="https://raw.githubusercontent.com/tobiasrausch/alfred/master/alfred.png">
   <h1></h1>
</p>

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/alfred/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/alfred/badges/downloads.svg)](https://anaconda.org/bioconda/alfred)
[![Build Status](https://travis-ci.org/tobiasrausch/alfred.svg?branch=master)](https://travis-ci.org/tobiasrausch/alfred)
[![Docker Build](https://img.shields.io/docker/build/trausch/alfred.svg)](https://hub.docker.com/r/trausch/alfred/)
[![GitHub license](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://raw.githubusercontent.com/tobiasrausch/alfred/master/LICENSE)
[![GitHub Releases](https://img.shields.io/github/release/tobiasrausch/alfred.svg)](https://github.com/tobiasrausch/alfred/releases)
[![GitHub Issues](https://img.shields.io/github/issues/tobiasrausch/alfred.svg)](https://github.com/tobiasrausch/alfred/issues)


Alfred installation
-------------------

Statically linked binaries are available from the [Alfred github release page](https://github.com/tobiasrausch/alfred/releases/). There is also an [Alfred Bioconda package](https://anaconda.org/bioconda/alfred).


Alfred installation from source
-------------------------------

To build Alfred from source you need some build essentials and the Boost libraries, i.e. for Ubuntu:

`apt install build-essential g++ cmake git-all liblzma-dev zlib1g-dev libbz2-dev liblzma-dev libboost-date-time-dev libboost-program-options-dev libboost-system-dev libboost-filesystem-dev libboost-iostreams-dev`

Once you have installed these system libraries you can compile and link Alfred.

`git clone --recursive https://github.com/tobiasrausch/alfred.git`

`cd alfred/`

`make all`

`make install`

`./bin/alfred -h`


BAM Alignment Quality Control
-----------------------------

Alfred computes various alignment metrics and summary statistics by read group

`./src/alfred qc -r <ref.fa> -o qc.tsv.gz <align.bam>`

Plotting alignment statistics

`Rscript scripts/stats.R qc.tsv.gz`

To convert all the alignment metrics from column format to row format for readability

`zgrep ^ME qc.tsv.gz | cut -f 2- | datamash transpose | column -t`


Interactive Quality Control Browser
-----------------------------------

Quality control metrics can be browsed interactively using the [web front end of Alfred](https://gear.embl.de/alfred).

`./src/alfred qc -r <ref.fa> -f json -o qc.json.gz <align.bam>`

Then just upload the qc.json.gz file to the Alfred GUI [https://gear.embl.de/alfred](https://gear.embl.de/alfred). A convenient feature of the web-front end is that multiple samples can be compared (as shown in the online example) if json files are merged prior to the upload.

`python ./scripts/merge.py sample1.json.gz sample2.json.gz sampleN.json.gz | gzip -c > multisample.json.gz`


BAM Alignment Quality Control for Targeted Sequencing
-----------------------------------------------------

If target regions are provided, Alfred computes the average coverage for each target and the on-target rate.

`./src/alfred qc -r <ref.fa> -b <targets.bed.gz> -o <qc.tsv.gz> <align.bam>`

For instance, for a human whole-exome data set.

`cd maps/ && Rscript exon.R`

`./src/alfred qc -r <hg19.fa> -b maps/exonic.hg19.bed.gz -o qc.tsv.gz <exome.bam>`

`Rscript scripts/stats.R qc.tsv.gz`

Alternatively, one can use the [interactive GUI](https://gear.embl.de/alfred) and upload the json file.

`./src/alfred qc -r <hg19.fa> -b maps/exonic.hg19.bed.gz -f json -o qc.json.gz <exome.bam>`


BAM Alignment Quality Control for ATAC-Seq
------------------------------------------

For ATAC-Seq data, the insert size distribution should reveal the DNA pitch and a clear nucleosome pattern with a peak for single nucleosomes and dimers. The transcription start site (TSS) enrichment should be >5 for a good ATAC-Seq library and ideally the duplicate rate is <20%, the alignment rate >70% and the standardized SD in coverage >0.3.

`cd maps/ && Rscript promoter.R`

`./src/alfred qc -r <hg19.fa> -b maps/hg19.promoter.bed.gz -o qc.tsv.gz <atac.bam>`

`Rscript scripts/stats.R qc.tsv.gz`

`zgrep ^ME qc.tsv.gz | datamash transpose | egrep "^Dup|^MappedFraction|^SD|^Enrich"`

ATAC-Seq libraries often have a large number of mitochondrial reads depending on the library preparation.

`zgrep ^CM qc.tsv.gz | egrep "Mapped|chrM"`

Alternatively, one can use the [interactive GUI](https://gear.embl.de/alfred) and upload the json file.

`./src/alfred qc -r <hg19.fa> -b maps/hg19.promoter.bed.gz -f json -o qc.json.gz <atac.bam>`


BAM Feature Counting for RNA-Seq
--------------------------------

Alfred can also assign reads to gene annotation features from a GTF file such as counting reads by gene or transcript identifier. Requires paired-end data.

`cd gtf/ && ./downloadGTF.sh`

`./src/alfred count_rna -g gtf/Homo_sapiens.GRCh37.75.gtf.gz <align.GRCh37.bam>`


BAM Read Counting for DNA-Seq
-----------------------------

For DNA sequencing, Alfred can be used to calculate the coverage in overlapping or non-overlapping windows or in given set of intervals. Requires paired-end data.

`./src/alfred count_dna -o <cov.gz> <align.GRCh37.bam>`

To plot the whole-chromosome coverage profile for chr1-22 and chrX.

`Rscript scripts/rd.R <cov.gz>`


BAM Feature Annotation
----------------------

Alfred can also be used to annotate peaks from ChIP-Seq or ATAC-Seq experiments. For instance to annotate overlapping/neighboring genes up to a distance of 10,000bp:

`./src/alfred annotate -d 10000 -g gtf/Homo_sapiens.GRCh37.75.gtf.gz <peaks.bed>`


Example E. coli data set
------------------------

The github source code includes a minimal example to check that alfred compiled properly from source and that the web front end is working.

`./src/alfred qc -r example/E.coli.fa.gz -o example/stats.tsv.gz example/E.coli.cram`

`Rscript scripts/stats.R example/stats.tsv.gz`

For the web front end.

`./src/alfred qc -r example/E.coli.fa.gz -f json -o ecoli.json.gz example/E.coli.cram`

Please upload ecoli.json.gz to the [Alfred web application](https://gear.embl.de/alfred).


Example plots
-------------

The [Alfred web application](https://gear.embl.de/alfred) includes a multi-sample example for whole-exome data including multiple read-groups. Additional examples for whole-genome sequencing, long-read sequencing, Hi-C, ATAC-Seq and other data sets are coming soon.


Credits
-------
[htseq-count](http://htseq.readthedocs.io) was used to benchmark and validate the RNA read counting.
