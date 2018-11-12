# Web application

Alfred's quality control JSON files can be interactively browsed with the
[companion web application](https://gear.embl.de/alfred).
All charts support panning and zooming and can be downloaded as PNG images.
The summary QC table can be downloaded as a CSV file.

To generate a quality control file in JSON format run [Alfred's command-line tool](/cli/) as follows:

```bash
alfred qc -r <ref.fa> -f json -o qc.json.gz <align.bam>
```

The output file `qc.json.gz` can then be uploaded at
[https://gear.embl.de/alfred/](https://gear.embl.de/alfred/).

## Features

An overview of all available charts and the most important alignment statistics provided by Alfred is below.

| Alignment Metric               | DNA-Seq (WGS) | DNA-Seq (Capture) | RNA-Seq | ChIP-Seq/ATAC-Seq | Chart Type         |
| ------------------------------ | ------------- | ----------------- | ------- | ----------------- | ------------------ |
| Mapping Statistics             | Yes           | Yes               | Yes     | Yes               | Table              |
| Duplicate Statistics           | Yes           | Yes               | Yes     | Yes               | Table              |
| Sequencing Error Rates         | Yes           | Yes               | Yes     | Yes               | Table              |
| Base Content Distribution      | Yes           | Yes               | Yes     | Yes               | Grouped Line Chart |
| Read Length Distribution       | Yes           | Yes               | Yes     | Yes               | Line Chart         |
| Base Quality Distribution      | Yes           | Yes               | Yes     | Yes               | Line Chart         |
| Coverage Histogram             | Yes           | Yes               | Yes     | Yes               | Line Chart         |
| Insert Size Distribution       | Yes           | Yes               | Yes     | Yes               | Grouped Line Chart |
| InDel Size Distribution        | Yes           | Yes               | Yes     | Yes               | Grouped Line Chart |
| InDel Context                  | Yes           | Yes               | Yes     | Yes               | Bar Chart          |
| GC Content                     | Yes           | Yes               | Yes     | Yes               | Grouped Line Chart |
| On-Target Rate                 |               | Yes               |         |                   | Line Chart         |
| Target Coverage Distribution   |               | Yes               |         |                   | Line Chart         |
| TSS Enrichment                 |               |                   |         | Yes               | Table              |
| DNA pitch / Nucleosome pattern |               |                   |         | Yes               | Grouped Line Chart |

## Base content distribution

The base content distribution shows any base calling bias along the read.
For an ideal library the lines for A, C, G, and T should run in parallel.
For a whole-genome assay the GC-content of that genome should be reflected in the relative amounts of each base.
Some libraries are expected to show a biased leading base distribution such as many RNA-Seq libraries
because of random hexamer priming or restriction based assays.

## Read length distribution

Illumina sequencers produce reads of fixed read length but long read technologies usually have a median read length >1000bp
and a long tail of reads with read lengths >30,000bp. This plot is also insightful to understand adapter trimming results
or the removal of low quality bases at the start or end of a read.

## Mean base quality distribution

This plot shows the mean base quality along the read. A typical Illumina profile shows base qualities >Q30
before base 30 and then a gradual loss of base quality accuracy towards the end of the read.

## Mapping quality distribution

This plot shows the mapping quality distribution for all mapped reads. The reported quality scores are aligner-dependent.

## Coverage histogram

The coverage histogram shows how many bases of the sequenced genome are at a given coverage.
Please note that for targeted assays (capture assays) this plot is expected to show a large portion of the genome at coverage=0.
For targeted assays, we therefore recommend checking the on-target rate and the targets above coverage level plots.

## On-target rate and targets above a given coverage level

For targeted assays, the two major concerns are capture efficiency (on-target rate)
and how many of the targets are ready for downstream analysis
(targets above a pre-defined coverage threshold).
A standard whole-exome sequencing assay yields at least 70% of reads on-target
(+/-200bp target extension) and at least 70% of targets >20x coverage.

## Insert size histogram

The insert size plot shows the outer insert size distribution for all read pairs
stratified by read pair orientation. There are different nomenclatures around for
defining the different paired-end layouts. The default Illumina paired-end layout is R+
(or forward-reverse, FR), the default Illumina mate-pair layout is R- (or reverse-forward, RF).
For specific sequencing assays, the insert size distribution can serve as a key quality control metric.
For instance, ATAC-Seq libraries should show the characteristic nucleosome pattern and DNA pitch.

## InDel size distribution

Histogram of indel sizes collected from all mapped reads. This plot aggregates the length
of all Cigar `I` and `D` operations.

## InDel Homopolymer Context

The homopolymer plot shows for all InDels (Cigar I and D operations) if the preceding 3 bases are
all A, all C, all G, or all T. If at least 2 different nucleotides are observed the reported
homopolymer context is "None". For Illumina reads, almost 50% of all reported InDels occur in a
homopolymer context with greater counts for A and T compared to G and C.

## GC content

To estimate a GC bias curve even for low-coverage single-cell data, Alfred computes for each mapped read
the local GC-content and then compares the estimated sample GC content to the expected, genome-wide GC content.
If a targeted assay is analyzed, Alfred, in addition, computes the GC content of all target regions.

## GC-Content and Mapping Statistics by Chromosome

This table lists the size, the number of Ns, the GC-content, and the number of mapped reads for each chromosome
as well as the observed-to-expected ratio of mapped reads.

## Summary statistics

The summary tab aggregates quality control data in a simple table that can be downloaded in CSV format.
This table is ideal to compare QC metrics across samples and/or sequencing assays. Among many other statistics,
the table lists, for instance, the number of duplicate reads, the number of unmapped reads, the number of
secondary and supplementary alignments, base-pair exact error rates stratified by mismatch,
insertion and deletion errors, and the median coverage and insert size of the sequenced sample.
The table provides more detailed statistics for specialized assays, i.e.
for 10X Genomics it lists the number of MI tagged reads, the total number of UMIs,
the fraction of haplotype-tagged reads and the N50 phased block length.
For ATAC-Seq data, users can provide a BED file of promoter regions and then the `EnrichmentOverBed` column
corresponds to TSS enrichment whereas for WES data, the enrichment quantifies the capturing efficiency
if the BED file contains all target regions.

## Example Data Sets

The [web application](https://gear.embl.de/alfred) hosts example data sets for a number of sequencing assays and sequencing technologies.

| Sequencing Assay  | Sequencing Technology |
| ----------------- | --------------------- |
| DNA-Seq (WGS)     | Illumina, PacBio, ONT |
| DNA-Seq (Capture) | Illumina              |
| RNA-Seq           | Illumina              |
| ATAC-Seq          | Illumina              |
| ChIP-Seq          | Illumina              |
