#!/bin/bash

if [ $# -lt 3 ] || [ $# -gt 4 ]; then
    echo "**********************************************************************"
    echo "Alfred"
    echo "This program comes with ABSOLUTELY NO WARRANTY."
    echo ""
    echo "Contact: Tobias Rausch (rausch@embl.de)"
    echo "**********************************************************************"
    echo ""
    echo "Usage: $0 <genome.fa.gz> <output prefix> <input.bam> [<targets.bed>]"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# Run analysis pipeline
if [ $# -eq 3 ]
then
    ${BASEDIR}/src/alfred -r ${1} -o ${2} ${3}
else
    ${BASEDIR}/src/alfred -r ${1} -b ${4} -o ${2} ${3}
fi

# Plot results
Rscript ${BASEDIR}/R/basequal.R ${2}.basequal.tsv
Rscript ${BASEDIR}/R/contentACGTN.R ${2}.contentACGTN.tsv
Rscript ${BASEDIR}/R/readlength.R ${2}.readlength.tsv
Rscript ${BASEDIR}/R/mapq.R ${2}.mapq.tsv
Rscript ${BASEDIR}/R/coverage.R ${2}.coverage.tsv
ICOL=`cat ${2}.metrics.tsv | head -n 1 | tr '\t' '\n' | awk '{print $0"\t"NR;}' | grep "^MedianInsertSize" | cut -f 2`
ISIZE=`cut -f 41 ${2}.metrics.tsv | tail -n +2 | awk '{SUM+=$1;} END {print SUM;}'`
FILES="${2}.basequal.tsv.png ${2}.contentACGTN.tsv.png ${2}.readlength.tsv.png ${2}.mapq.tsv.png"
if [ ${ISIZE} -ne 0 ]
then
    Rscript ${BASEDIR}/R/isize.R ${2}.isize.tsv
    FILES=${FILES}" ${2}.isize.tsv.png"
fi
if [ $# -eq 4 ]
then
    Rscript ${BASEDIR}/R/bedcov.R ${2}.bedcov.tsv
    Rscript ${BASEDIR}/R/ontarget.R ${2}.ontarget.tsv
    FILES=${FILES}" ${2}.bedcov.tsv.png ${2}.ontarget.tsv.png"
fi
Rscript R/metric.R ${2}.metrics.tsv
FILES=${FILES}" ${2}.metrics.tsv.pdf"

# Join plots
convert -adjoin ${FILES} ${2}.pdf
