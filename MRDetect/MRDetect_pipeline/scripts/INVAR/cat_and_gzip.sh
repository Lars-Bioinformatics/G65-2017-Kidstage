#!/bin/bash
#===================================================================================
#
# FILE: cat_and_gzip.sh
# USAGE: Do not run this script. Please refer to the Wiki.
#
# DESCRIPTION: Combine all final.ann.tsv files into one tsv.gz. 
#
#===================================================================================

# load config.sh
source config/config.sh

mkdir ${OUT}_gz

# awk
echo "checking that files with the right name are present"
ls ${OUT}/*BQ_${BASEQ}.MQ_${MAPQ}*final.ann.tsv

awk '$0 ~ /CHROM/ {next}; {print $0, FILENAME}' ${OUT}/*BQ_${BASEQ}.MQ_${MAPQ}*final.ann.tsv > ${OUT}_gz/${FINAL_PREFIX}.combined.final.ann.tsv

# replace tabs to preserve TSV structure
perl -p -i -e 's/ /\t/g' ${OUT}_gz/${FINAL_PREFIX}.combined.final.ann.tsv

echo gzipping
gzip -f ${OUT}_gz/${FINAL_PREFIX}.combined.final.ann.tsv
