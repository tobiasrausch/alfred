#!/bin/bash

if [ $# -ne 1 ]
then
    echo "Usage: $0 [tiny|full]"
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

if [ ${1} == "tiny" ]
then
    # Run the E.coli example
    echo "tiny example"
    ${BASEDIR}/../src/alfred qc -r ${BASEDIR}/../exampledata/E.coli.fa.gz -f json -o ecoli.json.gz ${BASEDIR}/../exampledata/E.coli.cram
elif [ ${1} == "full" ]
then
    echo "full example"

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
    python ${BASEDIR}/../scripts/merge.py > /dev/null
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

    if [ ! -f HG00733_I_045.alt_bwamem_GRCh38DH.20151004.CHS.monodinuc_strandseq.bam ]
    then
	wget 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20151203_strand_seq/HG00733/HG00733_I_045.alt_bwamem_GRCh38DH.20151004.CHS.monodinuc_strandseq.bam'
    fi
    if [ ! -f NA19238_I_026.alt_bwamem_GRCh38DH.20151004.CHS.monodinuc_strandseq.bam ]
    then
	wget 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20151203_strand_seq/NA19238/NA19238_I_026.alt_bwamem_GRCh38DH.20151004.CHS.monodinuc_strandseq.bam'
    fi
    ${BASEDIR}/../src/alfred qc -r GRCh38_full_analysis_set_plus_decoy_hla.fa -f json -o HG00733.json.gz HG00733_I_045.alt_bwamem_GRCh38DH.20151004.CHS.monodinuc_strandseq.bam
    ${BASEDIR}/../src/alfred qc -r GRCh38_full_analysis_set_plus_decoy_hla.fa -f json -o NA19238.json.gz NA19238_I_026.alt_bwamem_GRCh38DH.20151004.CHS.monodinuc_strandseq.bam
else
    echo "Unknown mode ${1}"
fi
