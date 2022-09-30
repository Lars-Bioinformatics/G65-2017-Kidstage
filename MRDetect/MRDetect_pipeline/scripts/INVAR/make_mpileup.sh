#!/bin/bash
#===================================================================================
#
# FILE: make_mpileup.sh
# USAGE: Do not run this file, please use INVAR1.sh. Please refer to the Wiki for guidance on running the pipeline.
#
# DESCRIPTION: The user should provide paths to MPILEUP, BCFTOOLS and FASTAREF if outside NRLAB.
#               These can be set in the config file
#
#
#===================================================================================

# read in variables from run_make_mpileup.sh
PREFIX=$1
SLOPPED_BED=$2
OUT=$3
NAME=$4
MAPQ=$5
BASEQ=$6

# load config.sh
source config/config.sh

TAGS="DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR"
[[ -z "$6" ]] && { echo "Please run run_make_mpileup.sh instead" ; exit 1; }

# samtools index each file
echo "running samtools index"
$SAMTOOLS_INDEX $PREFIX.bam

# remove .bam from input name
if [ -z ${dedup+x} ]; then
  echo "mpileup with --ff UNMAP"
  echo "using -R flag to ignore read groups"
  echo "${MPILEUP} -x -d ${MAX_DP} --ff UNMAP -q ${MAPQ} -Q ${BASEQ} -g -l ${SLOPPED_BED} -t ${TAGS} -f ${FASTAREF} ${PREFIX}.bam | ${BCFTOOLS} call -A -m -p 0 -o ${OUT}/${NAME}.vcf"
  ${MPILEUP} -x -R -d ${MAX_DP} --ff UNMAP -q ${MAPQ} -Q ${BASEQ} -g -l ${SLOPPED_BED} -t ${TAGS} -f ${FASTAREF} ${PREFIX}.bam | ${BCFTOOLS} call -A -m -p 0 -o ${OUT}/${NAME}.vcf
else
  echo "mpileup and deduplicate"
  echo "${MPILEUP} -x -d ${MAX_DP} -q ${MAPQ} -Q ${BASEQ} -g -l ${SLOPPED_BED} -t ${TAGS} -f ${FASTAREF} ${PREFIX}.bam | ${BCFTOOLS} call -A -m -p 0 -o ${OUT}/${NAME}.vcf"
  ${MPILEUP} -x -d ${MAX_DP} -q ${MAPQ} -Q ${BASEQ} -g -l ${SLOPPED_BED} -t ${TAGS} -f ${FASTAREF} ${PREFIX}.bam | ${BCFTOOLS} call -A -m -p 0 -o ${OUT}/${NAME}.vcf
fi


## split multiallelic sites into biallelic records
echo "splitting VCF"
${BCFTOOLS} norm -m -both ${OUT}/${NAME}.vcf -o ${OUT}/${NAME}.split.vcf

## get cols of interest
echo "obtaining cols of interest"
${BCFTOOLS} query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/DP4\t[%ADF]\t[%ADR]\t%INFO/MQSB\n' ${OUT}/${NAME}.split.vcf -o ${OUT}/${NAME}.split.tsv

awk -F"\t" 'BEGIN{OFS="\t"} { if($4==".") {print $1,$2,$3,$4,$5,$6,$7,"0",$8,"0",$9} else {split($7,a,","); split($8,b,","); print $1,$2,$3,$4,$5,$6,a[1],a[2],b[1],b[2],$9} }' ${OUT}/${NAME}.split.tsv > ${OUT}/${NAME}.final.tsv

## add slx+barcode as column
SLX=$(echo $NAME | grep -oP "SLX-[0-9]+")

if [ -z ${LIBRARY_PREP+x} ]; then
  echo "Don't recognize library prep will assume Rubicon"
  BC=$(echo $NAME |  grep -oP "D.*_D.*\\." | cut -d "." -f 1)

elif [ $LIBRARY_PREP == "Rubicon" ]; then
  BC=$(echo $NAME |  grep -oP "D.*_D.*\\." | cut -d "." -f 1)

elif [ $LIBRARY_PREP == "Agilent" ]; then
  BC=$(echo $NAME | grep -oP "SXTHS[A-Z]0[0-9]")

fi

sed -i "s/$/\t${SLX}/" ${OUT}/${NAME}.final.tsv
sed -i "s/$/\t${BC}/" ${OUT}/${NAME}.final.tsv

## add BQ and MQ information to the filename
echo "renaming file"
mv ${OUT}/${NAME}.final.tsv ${OUT}/${NAME}.BQ_$BASEQ.MQ_$MAPQ.final.tsv

#rm ${OUT}/${NAME}.vcf
#rm ${OUT}/${NAME}.split.vcf
#rm ${OUT}/${NAME}.split.tsv

# split tsv into 50 chunks
split -n l/${CHUNKS} -d ${OUT}/${NAME}.BQ_$BASEQ.MQ_$MAPQ.final.tsv temp/${NAME}.BQ_$BASEQ.MQ_$MAPQ.

## run annotation script
for FILE in $(ls temp/${NAME}.BQ_$BASEQ.MQ_$MAPQ.* | grep -v "ann" - ); do

  echo 'Processing:' $FILE

  sbatch --mem 10G \
    -o ${FILE}.ann.out \
    -e ${FILE}.ann.err \
    ${ANNOTATE_SCRIPT} $FILE;
done
