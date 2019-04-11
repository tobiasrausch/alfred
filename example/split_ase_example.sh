#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

echo "Download reference"
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa

echo "Download phased variants for chr22"
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi

echo "Subset to HG00732"
bcftools view -O b -o HG00732.bcf -s HG00732 -m2 -M2 -c 1 -C 1 ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz chr22:20000000-22000000
bcftools index HG00732.bcf

echo "Fetch HG00732 BAM slice"
samtools view -b -T GRCh38_full_analysis_set_plus_decoy_hla.fa ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/PUR/HG00732/alignment/HG00732.alt_bwamem_GRCh38DH.20150718.PUR.low_coverage.cram chr22:20000000-22000000 > HG00732.bam
samtools index HG00732.bam

echo "Split BAM"
${BASEDIR}/../src/alfred split -r GRCh38_full_analysis_set_plus_decoy_hla.fa -s HG00732 -v HG00732.bcf HG00732.bam

echo "Create allele-specific counts"
${BASEDIR}/../src/alfred ase -r GRCh38_full_analysis_set_plus_decoy_hla.fa -s HG00732 -v HG00732.bcf -p -f HG00732.bam
