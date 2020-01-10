You can build an [alfred](https://github.com/tobiasrausch/alfred) singularity container (SIF file) using

`sudo singularity build alfred.sif alfred.def`

Once you have built the container you can run analysis using

`singularity exec alfred.sif alfred qc -r ref.fa input.bam`
