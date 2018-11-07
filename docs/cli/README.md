# Usage

Alfred uses subcommands for [quality control](#alignment-quality-control) ([qc](#alignment-quality-control)), [feature counting](#bam-feature-counting) ([count_dna](#bam-read-counting-for-dna-seq), [count_rna](#bam-read-counting-for-rna-seq), [count_jct](#bam-read-counting-for-rna-seq)), [feature annotation](#bam-feature-annotation) ([annotate](#chip-seq-or-atac-seq-peak-annotation), tracks), alignment (pwalign, consensus) and haplotype-resolved analysis (split, ase). The subcommands are explained below.

## Alignment Quality Control

Alfred supports a command-line interface to run alignment quality control and a [web application](https://gear.embl.de) can be used to render all QC metrics.

### Command-line Interfact for BAM Quality Control

Alfred computes various alignment metrics and summary statistics by read group

```bash
alfred qc -r <ref.fa> -o qc.tsv.gz <align.bam>
```

Plotting alignment statistics

```bash
Rscript scripts/stats.R qc.tsv.gz
```

To convert all the alignment metrics from column format to row format for readability

```bash
zgrep ^ME qc.tsv.gz | cut -f 2- | datamash transpose | column -t
```

### Interactive Quality Control Browser

Quality control metrics can be browsed interactively using the [web front end of Alfred](https://gear.embl.de/alfred).

```bash
alfred qc -r <ref.fa> -f json -o qc.json.gz <align.bam>
```

Then just upload the qc.json.gz file to the Alfred GUI [https://gear.embl.de/alfred](https://gear.embl.de/alfred). A convenient feature of the web-front end is that multiple samples can be uploaded and compared.


### BAM Alignment Quality Control for Targeted Sequencing

If target regions are provided, Alfred computes the average coverage for each target and the on-target rate.

```bash
alfred qc -r <ref.fa> -b <targets.bed.gz> -o <qc.tsv.gz> <align.bam>
```

For instance, for a human whole-exome data set.

```bash
cd maps/ && Rscript exon.R
alfred qc -r <hg19.fa> -b maps/exonic.hg19.bed.gz -o qc.tsv.gz <exome.bam>
Rscript scripts/stats.R qc.tsv.gz
```

Alternatively, one can use the [interactive GUI](https://gear.embl.de/alfred) and upload the json file.

```bash
alfred qc -r <hg19.fa> -b maps/exonic.hg19.bed.gz -f json -o qc.json.gz <exome.bam>
```


### BAM Alignment Quality Control for ATAC-Seq

For ATAC-Seq data, the insert size distribution should reveal the DNA pitch and a clear nucleosome pattern with a peak for single nucleosomes and dimers. The transcription start site (TSS) enrichment should be >5 for a good ATAC-Seq library and ideally the duplicate rate is <20%, the alignment rate >70% and the standardized SD in coverage >0.3.

```bash
cd maps/ && Rscript promoter.R
alfred qc -r <hg19.fa> -b maps/hg19.promoter.bed.gz -o qc.tsv.gz <atac.bam>
Rscript scripts/stats.R qc.tsv.gz
zgrep ^ME qc.tsv.gz | datamash transpose | egrep "^Dup|^MappedFraction|^SD|^Enrich"
```

ATAC-Seq libraries often have a large number of mitochondrial reads depending on the library preparation.

```bash
zgrep ^CM qc.tsv.gz | egrep "Mapped|chrM"
```

Alternatively, one can use the [interactive GUI](https://gear.embl.de/alfred) and upload the json file.

```bash
alfred qc -r <hg19.fa> -b maps/hg19.promoter.bed.gz -f json -o qc.json.gz <atac.bam>
```

## BAM Feature Counting

Alfred supports counting reads in overlapping or non-overlapping windows, at predefined intervals in BED format,
or as gene and transcript counting for RNA-Seq in stranded or unstranded mode using a gtf or gff3 gene annotation
file. Expression values can be normalized as raw counts, FPKM, or FPKM-UQ values.

### BAM Read Counting for DNA-Seq

For DNA sequencing, Alfred can be used to calculate the coverage in overlapping or non-overlapping windows or in given set of intervals.

```bash
alfred count_dna -o <cov.gz> <align.GRCh37.bam>
```

To plot the whole-chromosome coverage profile for chr1-22 and chrX.

```bash
Rscript scripts/rd.R <cov.gz>
```

### BAM Read Counting for RNA-Seq

Alfred can also assign reads to gene annotation features from a GTF file such as counting reads by gene or transcript identifier.

```bash
cd gtf/ && ./downloadGTF.sh
alfred count_rna -g gtf/Homo_sapiens.GRCh37.75.gtf.gz <align.GRCh37.bam>
```

An experimental feature of Alfred is to count splice junction supporting reads. This method generates exon-exon junction counts for intra-gene exon-exon junctions, inter-gene exon-exon junctions and completely novel (not annotated) intra-chromosomal junctions.

```bash
alfred count_jct -g gtf/Homo_sapiens.GRCh37.75.gtf.gz <align.GRCh37.bam>
```

## BAM Feature Annotation

Alfred supports annotation of ChIP-Seq and ATAC-Seq peaks for neighboring genes or transcription factor binding sites (based on motif alignments). Additionally, browser tracks in UCSC bedgraph format can be computed with configurable resolution.

### ChIP-Seq or ATAC-Seq peak annotation

To annotate overlapping/neighboring genes up to a distance of 10,000bp:

```bash
alfred annotate -d 10000 -g gtf/Homo_sapiens.GRCh37.75.gtf.gz <peaks.bed>
```

The two output files summarize nearby genes by peak and vice versa (peaks by gene). Motif annotation based on alignment scores of motif-specific position weight matrices can also be obtained.

```bash
cd motif/ && ./downloadMotifs.sh
alfred annotate -r <hg19.fa> -m motif/jaspar.gz <peaks.bed>
```
