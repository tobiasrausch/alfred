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
	samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa
	cat GRCh38_full_analysis_set_plus_decoy_hla.fa | sed 's/v1//' | sed 's/v2//' > GRCh38_full_analysis_set_plus_decoy_hla.mod.fa
	samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.mod.fa
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
    if [ ! -f dna.exome.illumina.pe.ms.json.gz ]
    then
	for SAMPLE in HG00110 HG00111 HG00112 HG00113 HG00114 HG00115
	do
	    # Download 1000 Genomes exome cram file
	    if [ ! -f ${SAMPLE}.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram ]
	    then
		wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/GBR/${SAMPLE}/exome_alignment/${SAMPLE}.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram"
		wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/GBR/${SAMPLE}/exome_alignment/${SAMPLE}.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram.crai"
	    fi
	    ${BASEDIR}/../src/alfred qc -r GRCh38_full_analysis_set_plus_decoy_hla.fa -f json -o ${SAMPLE}.exome.illumina.pe.json.gz -b ${BASEDIR}/../maps/exonic.hg38.bed.gz ${SAMPLE}.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram
	done
	python ${BASEDIR}/../scripts/merge.py *.exome.illumina.pe.json.gz | gzip -c > dna.exome.illumina.pe.ms.json.gz
	rm *.exome.illumina.pe.json.gz
    fi

    # WGS
    if [ ! -f dna.wgs.illumina.pe.ms.json.gz ]
    then
	for SAMPLE in HG00512 HG00513
	do
	    if [ ! -f ${SAMPLE}.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.cram ]
	    then
		wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/data/CHS/${SAMPLE}/high_cov_alignment/${SAMPLE}.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.cram"
		wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/data/CHS/${SAMPLE}/high_cov_alignment/${SAMPLE}.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.cram.crai"
	    fi
	    ${BASEDIR}/../src/alfred qc -r GRCh38_full_analysis_set_plus_decoy_hla.fa -f json -o ${SAMPLE}.wgs.illumina.pe.json.gz ${SAMPLE}.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.cram
	done
	python ${BASEDIR}/../scripts/merge.py *.wgs.illumina.pe.json.gz | gzip -c > dna.wgs.illumina.pe.ms.json.gz
	rm *.wgs.illumina.pe.json.gz
    fi

    # WGS Mate-pairs
    if [ ! -f dna.wgs.illumina.mp.ms.json.gz ]
    then
	for SAMPLE in HG00512 HG00513
	do
	    if [ ! -f ${SAMPLE}.alt_bwamem_GRCh38DH.20150724.CHS.sv_7kb_mate.cram ]
	    then
		wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/data/CHS/${SAMPLE}/sv_7kb_mate/${SAMPLE}.alt_bwamem_GRCh38DH.20150724.CHS.sv_7kb_mate.cram"
		wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/data/CHS/${SAMPLE}/sv_7kb_mate/${SAMPLE}.alt_bwamem_GRCh38DH.20150724.CHS.sv_7kb_mate.cram.crai"
	    fi
	    ${BASEDIR}/../src/alfred qc -r GRCh38_full_analysis_set_plus_decoy_hla.fa -f json -o ${SAMPLE}.wgs.illumina.mp.json.gz ${SAMPLE}.alt_bwamem_GRCh38DH.20150724.CHS.sv_7kb_mate.cram
	done
	python ${BASEDIR}/../scripts/merge.py *.wgs.illumina.mp.json.gz | gzip -c > dna.wgs.illumina.mp.ms.json.gz
	rm *.wgs.illumina.mp.json.gz
    fi

    # RNA-Seq
    if [ ! -f rna.illumina.pe.ms.json.gz ]
    then
	for SAMPLE in SRS008746 SRS008747
	do
	    if [ ! -f ${SAMPLE}.bam ]
	    then
		wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20151026_strand_specific_mRNA/${SAMPLE}.gsnap_GRCh38Primary.20150922.PUR.mRNA.bam"
		wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20151026_strand_specific_mRNA/${SAMPLE}.gsnap_GRCh38Primary.20150922.PUR.mRNA.bam.bai"
		samtools view -H ${SAMPLE}.gsnap_GRCh38Primary.20150922.PUR.mRNA.bam | grep -v "ERCC" | grep -v "GL00" | grep -v "KI27" | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - ${SAMPLE}.gsnap_GRCh38Primary.20150922.PUR.mRNA.bam > ${SAMPLE}.bam
		samtools index ${SAMPLE}.bam
		rm ${SAMPLE}.gsnap_GRCh38Primary.20150922.PUR.mRNA.bam ${SAMPLE}.gsnap_GRCh38Primary.20150922.PUR.mRNA.bam.bai
	    fi
	    ${BASEDIR}/../src/alfred qc -r GRCh38_full_analysis_set_plus_decoy_hla.fa -f json -o ${SAMPLE}.rna.illumina.pe.json.gz ${SAMPLE}.bam
	done
	python ${BASEDIR}/../scripts/merge.py *.rna.illumina.pe.json.gz | gzip -c > rna.illumina.pe.ms.json.gz
	rm *.rna.illumina.pe.json.gz
    fi

    # PacBio
    if [ ! -f dna.pacbio.se.ms.json.gz ]
    then
	for SAMPLE in NA19239
	do
	    if [ ! -f ${SAMPLE}.pacbio.bam ]
	    then
		wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180102_pacbio_blasr_reheader/${SAMPLE}.pacbio-blasr-grch38-reheader.20180102.chr*.bam"
		for F in ${SAMPLE}.pacbio-blasr-grch38-reheader.20180102.*.bam
		do
		    if [ ! -f ${F}.bai ]
		    then
			samtools index ${F}
		    fi
		done
		samtools merge ${SAMPLE}.pacbio.bam ${SAMPLE}.pacbio-blasr-grch38-reheader.20180102.*.bam
		samtools index ${SAMPLE}.pacbio.bam
		rm ${SAMPLE}.pacbio-blasr-grch38-reheader.20180102.chr*
	    fi
	    ${BASEDIR}/../src/alfred qc -r GRCh38_full_analysis_set_plus_decoy_hla.mod.fa -f json -o ${SAMPLE}.dna.pacbio.se.json.gz ${SAMPLE}.pacbio.bam
	done
	mv ${SAMPLE}.dna.pacbio.se.json.gz dna.pacbio.se.ms.json.gz 
    fi
    
    # Hi-C
    if [ ! -f hic.illumina.pe.ms.json.gz ]
    then
	for SAMPLE in HG00732 HG00733
	do
	    if [ ! -f ${SAMPLE}_Hi-C_biorep2_merged_filtered.bam ]
	    then
		wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20160822_HiC_bam_files/${SAMPLE}_Hi-C_biorep2_merged_filtered.bam"
		samtools index ${SAMPLE}_Hi-C_biorep2_merged_filtered.bam
	    fi
	    ${BASEDIR}/../src/alfred qc -r GRCh38_full_analysis_set_plus_decoy_hla.fa -f json -o ${SAMPLE}.hic.illumina.pe.json.gz ${SAMPLE}_Hi-C_biorep2_merged_filtered.bam
	done
	python ${BASEDIR}/../scripts/merge.py *.hic.illumina.pe.json.gz | gzip -c > hic.illumina.pe.ms.json.gz
	rm *.hic.illumina.pe.json.gz
    fi

    # ONT
    if [ ! -f dna.wgs.ont.se.ms.json.gz ]
    then
	if [ -f /opt/dev/HG00733/HG00733.bam ]
	then
	    ${BASEDIR}/../src/alfred qc -r GRCh38_full_analysis_set_plus_decoy_hla.fa -f json -o dna.wgs.ont.se.ms.json.gz /opt/dev/HG00733/HG00733.bam
	fi
    fi
else
    echo "Unknown mode ${1}"
fi

