# Web application

Alfred's quality control json files can be interactively browsed with the companion [web application](https://gear.embl.de/alfred). All charts support panning and zooming. Charts can be downloaded as png. Tables of quality control metrics can be downloaded as csv. To generate a quality control file in json format:

```bash
alfred qc -r <ref.fa> -f json -o qc.json.gz <align.bam>
```

`qc.json.gz` can then be uploaded at [https://gear.embl.de/alfred](https://gear.embl.de/alfred).



## Features

An overview of all available charts and the most important alignment statistics provided by Alfred is below.


| Alignment Metric               | DNA-Seq (WGS)  | DNA-Seq (Capture)  | RNA-Seq  |  ChIP-Seq/ATAC-Seq | Chart Type |
|--------------------------------|----------------|--------------------|----------|--------------------|------------|
| Mapping Statistics             | Yes | Yes  | Yes | Yes  | Table|
| Duplicate Statistics           | Yes | Yes  | Yes | Yes  | Table|
| Sequencing Error Rates         | Yes | Yes  | Yes | Yes  | Table|
| Base Content Distribution      | Yes | Yes  | Yes | Yes  | Grouped Line Chart|
| Read Length Distribution       | Yes | Yes  | Yes | Yes  | Line Chart|
| Base Quality Distribution      | Yes | Yes  | Yes | Yes  | Line Chart|
| Coverage Histogram             | Yes | Yes  | Yes | Yes  | Line Chart|
| Insert Size Distribution       | Yes | Yes  | Yes | Yes  | Grouped Line Chart|
| InDel Size Distribution        | Yes | Yes  | Yes | Yes  | Grouped Line Chart|
| InDel Context                  | Yes | Yes  | Yes | Yes  | Bar Chart|
| GC Content                     | Yes | Yes  | Yes | Yes  | Grouped Line Chart|
| On-Target Rate                 |     | Yes  |     |      | Line Chart|
| Target Coverage Distribution   |     | Yes  |     |      | Line Chart|
| TSS Enrichment                 |     |      |     | Yes  | Table|
| DNA pitch / Nucleosome pattern |     |      |     | Yes  | Grouped Line Chart|


## Example Data Sets

The [web application](https://gear.embl.de/alfred) hosts example data sets for a number of sequencing assays and sequencing technologies.

| Sequencing Assay   | Sequencing Technology |
|--------------------|-----------------------|
| DNA-Seq (WGS)      | Illumina, PacBio, ONT |
| DNA-Seq (Capture)  | Illumina              |
| RNA-Seq            | Illumina              |
| ATAC-Seq           | Illumina              |
| ChIP-Seq           | Illumina              |
