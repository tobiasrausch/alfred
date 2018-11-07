# FAQ

::: tip
For questions, help or feature requests please contact gear_genomics@embl.de
:::

[[toc]]

## Does Alfred support the CRAM format?

Yes, Alfred uses [HTSlib](https://github.com/samtools/htslib) to read/write BAM/CRAM files.

## Is there an example data set to test my Alfred installation?

The github source code includes a minimal example to check that alfred compiled properly from source and that the web front end is working.

```bash
alfred qc -r example/E.coli.fa.gz -o example/stats.tsv.gz example/E.coli.cram
Rscript scripts/stats.R example/stats.tsv.gz
```

For the web front end.

```bash
alfred qc -r example/E.coli.fa.gz -f json -o ecoli.json.gz example/E.coli.cram
```

Please upload `ecoli.json.gz` to the [Alfred web application](https://gear.embl.de/alfred).

## Is the feature counting paired-end aware?

Yes, Alfred counts fragments (read pairs) and not individual reads.

## Why are hard clipping statistics always zero?

Many aligners trim primary alignments using soft-clips and only secondary and supplementary alignments use hard clips. For long reads you may want to evaluate secondary and supplementary alignments using the `-s` and `-u` command-line flags.

```bash
alfred qc -su -r <genome.fa> <input.bam>
```
