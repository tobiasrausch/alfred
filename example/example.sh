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
    ${BASEDIR}/../src/alfred qc -r ${BASEDIR}/E.coli.fa.gz -f json -o ecoli.json.gz ${BASEDIR}/E.coli.cram
    ${BASEDIR}/../src/alfred qc -r ${BASEDIR}/E.coli.fa.gz -o ecoli.tsv.gz ${BASEDIR}/E.coli.cram
    Rscript ${BASEDIR}/../scripts/stats.R ecoli.tsv.gz

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
    if [ ! -f ${BASEDIR}/../maps/exonic.hg38.bed.gz ]
    then
	cd ${BASEDIR}/../maps/ && Rscript exon.R && cd ${BASEDIR}
    fi
    
    # Download 1000 Genomes exome cram file
    if [ ! -f HG00114.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram ]
    then
	wget 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/GBR/HG00114/exome_alignment/HG00114.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram'
	wget 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/GBR/HG00114/exome_alignment/HG00114.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram.crai'
    fi

    # Run alfred
    ${BASEDIR}/../src/alfred qc -r GRCh38_full_analysis_set_plus_decoy_hla.fa -f json -o HG00114.exome.illumina.json.gz -b ${BASEDIR}/../maps/exonic.hg38.bed.gz HG00114.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram
    ${BASEDIR}/../src/alfred qc -r GRCh38_full_analysis_set_plus_decoy_hla.fa -o HG00114.exome.illumina.tsv.gz -b ${BASEDIR}/../maps/exonic.hg38.bed.gz HG00114.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram
    Rscript ${BASEDIR}/../scripts/stats.R HG00114.exome.illumina.tsv.gz

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
	if [ ! -f ${SAMPLE}.exome.illumina.json.gz ]
	then
	    ${BASEDIR}/../src/alfred qc -r GRCh38_full_analysis_set_plus_decoy_hla.fa -f json -o ${SAMPLE}.exome.illumina.json.gz -b ${BASEDIR}/../maps/exonic.hg38.bed.gz ${SAMPLE}.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram
	fi
    done
    python ${BASEDIR}/../scripts/merge.py *.exome.illumina.json.gz | gzip -c > dna.exome.illumina.ms.json.gz

    # Low coverage WGS
    for SAMPLE in HG00110 HG00111 HG00112 HG00113 HG00114 HG00115
    do
	if [ ! -f ${SAMPLE}.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage.cram ]
	then
	    wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/GBR/${SAMPLE}/alignment/${SAMPLE}.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage.cram"
	    wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/GBR/${SAMPLE}/alignment/${SAMPLE}.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage.cram.crai"
	fi
	if [ ! -f ${SAMPLE}.wgs.illumina.json.gz ]
	then
	    ${BASEDIR}/../src/alfred qc -r GRCh38_full_analysis_set_plus_decoy_hla.fa -f json -o ${SAMPLE}.wgs.illumina.json.gz ${SAMPLE}.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage.cram
	fi
    done
    python ${BASEDIR}/../scripts/merge.py *.wgs.illumina.json.gz | gzip -c > dna.wgs.illumina.ms.json.gz

    # PacBio
    for SAMPLE in NA19239 NA19238
    do
	if [ ! -f ${SAMPLE}.pacbio.bam ]
	then
	    wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180102_pacbio_blasr_reheader/${SAMPLE}.pacbio-blasr-grch38-reheader.20180102.chr*.bam.bai"
	    wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180102_pacbio_blasr_reheader/${SAMPLE}.pacbio-blasr-grch38-reheader.20180102.chr*.bam"
	    samtools merge ${SAMPLE}.pacbio.bam ${SAMPLE}.pacbio-blasr-grch38-reheader.20180102.*.bam
	    samtools index ${SAMPLE}.pacbio.bam
	fi
    done
    
    # Hi-C
    for SAMPLE in HG00732 HG00733
    do
	if [ ! -f ${SAMPLE}_Hi-C_biorep2_merged_filtered.bam ]
	then
	    wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20160822_HiC_bam_files/${SAMPLE}_Hi-C_biorep2_merged_filtered.bam"
	    samtools index ${SAMPLE}_Hi-C_biorep2_merged_filtered.bam
	fi
    done
else
    echo "Unknown mode ${1}"
fi

