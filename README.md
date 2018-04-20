<p align="center">
<img width="230" height="100" src="https://raw.githubusercontent.com/tobiasrausch/alfred/master/alfred.png">
</p>

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/alfred/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/alfred/badges/downloads.svg)](https://anaconda.org/bioconda/alfred)
[![Build Status](https://travis-ci.org/dellytools/delly.svg?branch=master)](https://travis-ci.org/tobiasrausch/alfred)
[![Docker Automated buil](https://img.shields.io/docker/automated/jrottenberg/ffmpeg.svg?style=flat-square)](https://hub.docker.com/r/trausch/alfred/)


Alfred installation
---------------------

The easiest way to get Alfred is to download a statically linked binary from the [Alfred github release page](https://github.com/tobiasrausch/alfred/releases/).
Alternatively, you can build Alfred from source. Alfred dependencies are included as submodules so you need to do a recursive clone.

`git clone --recursive https://github.com/tobiasrausch/alfred.git`

`cd alfred/`

`make all`

BAM Alignment Quality Control
-----------------------------

Alfred computes various alignment metrics and summary statistics by read group

`./src/alfred qc -r <ref.fa> -o <outfile.tsv.gz> <align.bam>`

Plotting alignment statistics

`Rscript R/stats.R <outfile.tsv.gz>`

To convert all the alignment metrics from column format to rows to easily read it on screen

`zgrep ^ME stats.tsv.gz | cut -f 2- | datamash transpose | column -t`


BAM Alignment Quality Control for Targeted Sequencing
-----------------------------------------------------

If target regions are provided, Alfred computes the average coverage for each target and the on-target rate.

`./src/alfred qc -r <ref.fa> -b <targets.bed> -o <outfile.tsv.gz> <align.bam>`

For instance, for a human whole-exome data set.

`cd maps/ && Rscript exon.R`

`./src/alfred qc -r <hg19.fa> -b maps/exonic.hg19.bed.gz -o <outfile.tsv.gz> <exome.bam>`

`Rscript R/stats.R <outfile.tsv.gz>`


BAM Alignment Quality Control for ATAC-Seq
------------------------------------------

For ATAC-Seq data, the insert size distribution should reveal the DNA pitch and a clear nucleosome pattern with a peak for single nucleosomes and dimers. The transcription start site (TSS) enrichment should be >5 for a good ATAC-Seq library and ideally the duplicate rate is <20%, the alignment rate >70% and the standardized SD in coverage >0.3.

`cd maps/ && Rscript promoter.R`

`./src/alfred qc -r <hg19.fa> -b maps/hg19.promoter.bed.gz -o <outfile.tsv.gz> <atac.bam>`

`Rscript R/stats.R <outfile.tsv.gz>`

`zgrep ^ME <outfile.tsv.gz> | datamash transpose | egrep "DuplicateFraction|MappedFraction|SDCoverage|EnrichmentOverBed"`

ATAC-Seq often has a high number of mitochondrial reads depending on the library preparation.




Web Front End
-------------

The Genome Analysis Server [GEAR](https://gear.embl.de) has a [web front end](https://gear.embl.de/alfred) for Alfred. Instead of plotting the alignment metrics and statistics with R you can also upload the QC output file of Alfred [here](https://gear.embl.de/alfred).



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

`Rscript R/rd.R <cov.gz>`


BAM Feature Annotation
----------------------

Alfred can also be used to annotate peaks from ChIP-Seq or ATAC-Seq experiments. For instance to annotate overlapping/neighboring genes up to a distance of 10,000bp:

`./src/alfred annotate -d 10000 -g gtf/Homo_sapiens.GRCh37.75.gtf.gz <peaks.bed>`


Example E. coli data set
------------------------

`./src/alfred qc -r exampledata/E.coli.fa.gz -o exampledata/stats.tsv.gz exampledata/E.coli.cram`

`Rscript R/stats.R exampledata/stats.tsv.gz`

The final pdf is exampledata/stats.tsv.gz.pdf


Example QC plots
----------------

[Whole-genome paired-end data with multiple read groups](https://raw.githubusercontent.com/tobiasrausch/alfred/master/exampleplots/NA06985.pe.pdf)

[Whole-exome paired-end data with on-target rate](https://raw.githubusercontent.com/tobiasrausch/alfred/master/exampleplots/HG00112.wes.pdf)

[Whole-genome jumping library](https://raw.githubusercontent.com/tobiasrausch/alfred/master/exampleplots/HG00513.mp.pdf)


Credits
-------
[htseq-count](http://htseq.readthedocs.io) was used to benchmark and validate the RNA read counting.
