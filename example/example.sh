#!/bin/bash

if [ $# -ne 1 ]
then
    echo "Usage: $0 [tiny|full]"
    exit -1
fi

# Check dependency tools
wget --version > /dev/null
if [ $? -ne 0 ]
then
    echo ""
    echo "wget is required!"
    echo ""
    exit
fi
samtools --version > /dev/null
if [ $? -ne 0 ]
then
    echo ""
    echo "Samtools is required!"
    echo ""
    exit
fi
python ../scripts/merge.py > /dev/null
if [ $? -ne 0 ]
then
    echo ""
    echo "Python is required!"
    echo ""
    exit
fi

# Download and index reference
if [ ! -f GRCh38_full_analysis_set_plus_decoy_hla.fa ]
then
    wget 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa'
fi
if [ ! -f GRCh38_full_analysis_set_plus_decoy_hla.fa.fai ]
then
    samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa
fi

# Download data
if [ ${1} == "full" ]
then
    echo "full example"
else
    echo "tiny example"
    if [ ! -f HG00733_I_045.alt_bwamem_GRCh38DH.20151004.CHS.monodinuc_strandseq.bam ]
    then
	wget 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20151203_strand_seq/HG00733/HG00733_I_045.alt_bwamem_GRCh38DH.20151004.CHS.monodinuc_strandseq.bam'
    fi
    if [ ! -f NA19238_I_026.alt_bwamem_GRCh38DH.20151004.CHS.monodinuc_strandseq.bam ]
    then
	wget 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20151203_strand_seq/NA19238/NA19238_I_026.alt_bwamem_GRCh38DH.20151004.CHS.monodinuc_strandseq.bam'
    fi
    ../src/alfred qc -r GRCh38_full_analysis_set_plus_decoy_hla.fa -f json -o HG00733.json.gz HG00733_I_045.alt_bwamem_GRCh38DH.20151004.CHS.monodinuc_strandseq.bam
    ../src/alfred qc -r GRCh38_full_analysis_set_plus_decoy_hla.fa -f json -o NA19238.json.gz NA19238_I_026.alt_bwamem_GRCh38DH.20151004.CHS.monodinuc_strandseq.bam
fi
