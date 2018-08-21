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
    
    # Download and index reference
    if [ ! -f GRCh38_full_analysis_set_plus_decoy_hla.fa ]
    then
	wget 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa'
    fi
    if [ ! -f GRCh38_full_analysis_set_plus_decoy_hla.fa.fai ]
    then
	samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa
    fi

    # Download gene annotation
    if [ ! -f ${BASEDIR}/../gtf/Homo_sapiens.GRCh38.91.gtf.gz ]
    then
	cd ${BASEDIR}/../gtf/ && ./downloadGTF.sh && cd ${BASEDIR}
    fi

    # Generate exon target file
    if [ ! -f hg38.exon.bed.gz ]
    then
 	zcat ${BASEDIR}/../gtf/Homo_sapiens.GRCh38.91.gtf.gz  | grep "^[1-9X]" | grep -P "\texon\t" | grep "protein_coding" | grep -P "\tensembl" | cut -f 1,4,5 | sort -k1,1V -k2,2n | sed 's/^/chr/' | gzip -c > hg38.exon.bed.gz
    fi
    
    # Download 1000 Genomes exome cram file
    if [ ! -f HG00114.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram ]
    then
	wget 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/GBR/HG00114/exome_alignment/HG00114.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram'
	wget 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/GBR/HG00114/exome_alignment/HG00114.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram.crai'
    fi

    # Run alfred
    ${BASEDIR}/../src/alfred qc -r GRCh38_full_analysis_set_plus_decoy_hla.fa -f json -o HG00114.exome.json.gz -b hg38.exon.bed.gz HG00114.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram

elif [ ${1} == "benchmark" ]
then
    ########
    # Requires a "full" run first to download required annotation and reference files
    ########

    # Exome
    for SAMPLE in HG00110 HG00111 HG00112 HG00113 HG00114 HG00115
    do
	# Download 1000 Genomes exome cram file
	if [ ! -f ${SAMPLE}.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram ]
	then
	    wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/GBR/${SAMPLE}/exome_alignment/${SAMPLE}.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram"
	    wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/GBR/${SAMPLE}/exome_alignment/${SAMPLE}.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram.crai"
	fi
	${BASEDIR}/../src/alfred qc -r GRCh38_full_analysis_set_plus_decoy_hla.fa -f json -o ${SAMPLE}.exome.json.gz -b hg38.exon.bed.gz ${SAMPLE}.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram
    done

    # Low coverage WGS
    for SAMPLE in HG00110 HG00111 HG00112 HG00113 HG00114 HG00115
    do
	if [ ! -f ${SAMPLE}.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage.cram ]
	then
	    wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/GBR/${SAMPLE}/alignment/${SAMPLE}.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage.cram"
	    wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/GBR/${SAMPLE}/alignment/${SAMPLE}.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage.cram.crai"
	fi
	${BASEDIR}/../src/alfred qc -r GRCh38_full_analysis_set_plus_decoy_hla.fa -f json -o ${SAMPLE}.wgs.json.gz ${SAMPLE}.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage.cram
    done

else
    echo "Unknown mode ${1}"
fi

