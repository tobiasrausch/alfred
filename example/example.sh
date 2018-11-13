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
    ${BASEDIR}/../src/alfred qc -r ${BASEDIR}/E.coli.fa.gz -j ecoli.json.gz ${BASEDIR}/E.coli.cram
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
	wget 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz'
	zcat human_g1k_v37.fasta.gz  | sed 's/>\([0-9XYM][0-9T]*\) />chr\1 /' | sed 's/>chrMT/>chrM/' > hg19.fa
	rm human_g1k_v37.fasta.gz
	samtools faidx hg19.fa
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
    ${BASEDIR}/../src/alfred qc -r GRCh38_full_analysis_set_plus_decoy_hla.fa -j HG00114.exome.illumina.json.gz -b ${BASEDIR}/../maps/exonic.hg38.bed.gz HG00114.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram
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
	    ${BASEDIR}/../src/alfred qc -r GRCh38_full_analysis_set_plus_decoy_hla.fa -j ${SAMPLE}.exome.illumina.pe.json.gz -b ${BASEDIR}/../maps/exonic.hg38.bed.gz ${SAMPLE}.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram
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
	    ${BASEDIR}/../src/alfred qc -r GRCh38_full_analysis_set_plus_decoy_hla.fa -j ${SAMPLE}.wgs.illumina.pe.json.gz ${SAMPLE}.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.cram
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
	    ${BASEDIR}/../src/alfred qc -r GRCh38_full_analysis_set_plus_decoy_hla.fa -j ${SAMPLE}.wgs.illumina.mp.json.gz ${SAMPLE}.alt_bwamem_GRCh38DH.20150724.CHS.sv_7kb_mate.cram
	done
	python ${BASEDIR}/../scripts/merge.py *.wgs.illumina.mp.json.gz | gzip -c > dna.wgs.illumina.mp.ms.json.gz
	rm *.wgs.illumina.mp.json.gz
    fi

    # RNA-Seq, Geuvadis
    if [ ! -f rna.illumina.pe.ms.json.gz ]
    then
	for SAMPLE in HG00096.1.M_111124_6 HG00101.1.M_111124_4 HG00104.1.M_111124_5 HG00117.1.M_111124_2 HG00121.1.M_111124_7 
	do
	    if [ ! -f ${SAMPLE}.bam ]
	    then
		wget "https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/${SAMPLE}.bam"
		samtools index ${SAMPLE}.bam
	    fi
	    ${BASEDIR}/../src/alfred qc -a ${SAMPLE} -r hg19.fa -j ${SAMPLE}.rna.illumina.pe.json.gz ${SAMPLE}.bam
	done
	python ${BASEDIR}/../scripts/merge.py *.rna.illumina.pe.json.gz | gzip -c > rna.illumina.pe.ms.json.gz
	rm *.rna.illumina.pe.json.gz
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
	    ${BASEDIR}/../src/alfred qc -r GRCh38_full_analysis_set_plus_decoy_hla.fa -j ${SAMPLE}.hic.illumina.pe.json.gz ${SAMPLE}_Hi-C_biorep2_merged_filtered.bam
	done
	python ${BASEDIR}/../scripts/merge.py *.hic.illumina.pe.json.gz | gzip -c > hic.illumina.pe.ms.json.gz
	rm *.hic.illumina.pe.json.gz
    fi

    # ONT
    if [ ! -f dna.wgs.ont.se.ms.json.gz ]
    then
	if [ -f /opt/dev/HG00733/HG00733.bam ]
	then
	    ${BASEDIR}/../src/alfred qc -r GRCh38_full_analysis_set_plus_decoy_hla.fa -j dna.wgs.ont.se.ms.json.gz /opt/dev/HG00733/HG00733.bam
	fi
    fi

    # PacBio
    if [ ! -f dna.wgs.pacbio.se.ms.json.gz ]
    then
	for SAMPLE in NA19238 NA19239
	do
	    if [ ! -f ${SAMPLE}_bwamem_GRCh38DH_YRI_20160905_pacbio.bam ]
	    then
		wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20160905_smithm_pacbio_aligns/${SAMPLE}_bwamem_GRCh38DH_YRI_20160905_pacbio.bam"
		wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20160905_smithm_pacbio_aligns/${SAMPLE}_bwamem_GRCh38DH_YRI_20160905_pacbio.bam.bai"
	    fi
	    ${BASEDIR}/../src/alfred qc -r GRCh38_full_analysis_set_plus_decoy_hla.fa -j ${SAMPLE}.dna.wgs.pacbio.se.json.gz ${SAMPLE}_bwamem_GRCh38DH_YRI_20160905_pacbio.bam
	done
	python ${BASEDIR}/../scripts/merge.py *.dna.wgs.pacbio.se.json.gz | gzip -c > dna.wgs.pacbio.se.ms.json.gz
	rm *.dna.wgs.pacbio.se.json.gz
    fi

    # ATAC-Seq
    if [ ! -f atac.illumina.pe.ms.json.gz ]
    then
	for SAMPLE in atac1 atac2
	do
	    ${BASEDIR}/../src/alfred qc -a ${SAMPLE} -r hs37d5.fa -j ${SAMPLE}.atac.illumina.pe.json.gz ${SAMPLE}.bam
	done
	python ${BASEDIR}/../scripts/merge.py *.atac.illumina.pe.json.gz | gzip -c > atac.illumina.pe.ms.json.gz
	rm *.atac.illumina.pe.json.gz
    fi

    # ChIP-Seq, Encode
    if [ ! -f chip.illumina.se.ms.json.gz ]
    then
	for SAMPLE in H3K27ac  H3K27me3  H3K9ac
	do
	    ${BASEDIR}/../src/alfred qc -a GM12878_${SAMPLE} -r hs37d5.fa -j ${SAMPLE}.chip.illumina.se.json.gz ${SAMPLE}.bam
	done
	python ${BASEDIR}/../scripts/merge.py *.chip.illumina.se.json.gz | gzip -c > chip.illumina.se.ms.json.gz
	rm *.chip.illumina.se.json.gz
    fi
elif [ ${1} == "runtime" ]
then
    # RNA benchmark
    zcat ${BASEDIR}/../gtf/Homo_sapiens.GRCh37.75.gtf.gz | grep "^#" | gzip -c > Homo_sapiens.GRCh37.75.chr.gtf.gz
    zcat ${BASEDIR}/../gtf/Homo_sapiens.GRCh37.75.gtf.gz  | grep -v "^#" | grep -P "^[0-9XY]*\t" | sed 's/^/chr/' | gzip -c >> Homo_sapiens.GRCh37.75.chr.gtf.gz
    /usr/bin/time -v ${BASEDIR}/../src/alfred count_rna -g Homo_sapiens.GRCh37.75.chr.gtf.gz HG00104.1.M_111124_5.bam
    gunzip Homo_sapiens.GRCh37.75.chr.gtf.gz
    /usr/bin/time -v ./qualimap_v2.2.1/qualimap comp-counts -pe -id gene_id -type exon -out qualimap.count -gtf Homo_sapiens.GRCh37.75.chr.gtf -bam HG00104.1.M_111124_5.bam
    /usr/bin/time -v htseq-count -s no -f bam -r pos HG00104.1.M_111124_5.bam Homo_sapiens.GRCh37.75.chr.gtf > htseq.count
    rm Homo_sapiens.GRCh37.75.chr.gtf

    # DNA benchmark
    zcat ${BASEDIR}/../maps/exonic.hg19.bed.gz | sed 's/^chr//' | gzip -c > exonic.hg19.bed.gz
    # QC
    /usr/bin/time -v ${BASEDIR}/../src/alfred qc -j qc.json.gz -b exonic.hg19.bed.gz -r 1kGP.fa HG00111.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam
    /usr/bin/time -v ./qualimap_v2.2.1/qualimap bamqc --java-mem-size=4G -nt 1 -bam HG00111.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam -sd -gd "HUMAN - hg19" -outformat PDF:HTML
    # Window counting
    /usr/bin/time -v ${BASEDIR}/../src/alfred count_dna -i exonic.hg19.bed.gz -m 0 HG00111.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam
    /usr/bin/time -v bedtools multicov -p -bams HG00111.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam -bed exonic.hg19.bed.gz  > bedtools.count
    rm exonic.hg19.bed.gz

    # Browser Tracks
    /usr/bin/time -v ${BASEDIR}/../src/alfred tracks HG00111.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam
    /usr/bin/time -v makeTagDirectory tagdir -genome 1kGP.fa HG00111.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam
    /usr/bin/time -v makeUCSCfile tagdir -style dnase -fsize 5e7 -o homer.bedGraph
else
    echo "Unknown mode ${1}"
fi

