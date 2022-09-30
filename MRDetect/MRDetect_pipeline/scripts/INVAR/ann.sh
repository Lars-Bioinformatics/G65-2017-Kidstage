#!/bin/bash
#===================================================================================
#
# FILE: ann.sh
# USAGE: do not run this .sh. Please refer to the Wiki for guidance on running the code.
#
# DESCRIPTION: This shell script contains code to annotate mpileup data with trinucleotide context
#               and COSMIC mutation information
#
#===================================================================================
FILE=$1

# load config.sh
source config/config.sh

echo "running annotation"

# remove loci with MQSB == "."
echo "awk"
awk -v MIN_DP=$MIN_DP '$5 >= MIN_DP {print $0}' $FILE > $FILE.awk.tsv

sed 1i"CHROM\tPOS\tREF\tALT\tDP\tDP4\tREF_F\tALT_F\tREF_R\tALT_R\tMQSB\tSLX\tBARCODE" $FILE.awk.tsv > $FILE.header.tsv

# run trinucleotide context and annotate with SNPs and Cosmic
echo "running trinucleotide context"
  $ANNOTATE_PYTHON \
  $K1G_DB \
  $COSMIC_DB $FILE.header.tsv > $FILE.ann.tsv

# remove temp files
echo "removing temp files"
rm $FILE.header.tsv
rm $FILE.awk.tsv

# moving output to output folder
echo "moving files"
FILE_PREFIX=$(echo $FILE | cut -d'/' -f 2)
echo "FILE_PREFIX = ${FILE_PREFIX}"
mv $FILE.ann.tsv ${OUT}/${FILE_PREFIX}.final.ann.tsv
