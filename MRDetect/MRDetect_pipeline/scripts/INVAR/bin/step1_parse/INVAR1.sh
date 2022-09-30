#!/bin/bash
#===================================================================================
#
# FILE: INVAR1.sh
# USAGE: ./INVAR1.sh
#
# DESCRIPTION: Run mpileup on error-suppressed BAM files from individualised capture data.
#              mpileup data can then be read into R for INVAR analysis.
#              The user needs to define paths to software (mpileup, bcftools)
#              in the make_mpileup.sh script.
#
# REQUIREMENTS: samtools 1.3.1, bcftools 1.3.1
#
#===================================================================================

rm -rf temp
mkdir temp
#rm log.txt

chmod -R 775 ~/tapas_pipeline

# load config.sh
echo "loading config.sh"
source config/config.sh

if [[ ! -f $INPUT_FILES ]] ; then
    echo 'Input file txt does not exist, aborting.'
    echo 'This txt file contains a list of paths to BAM files to be run.'
    exit
fi


if [ ! -d "$OUT" ]; then
mkdir -p $OUT
fi

regex="^(.+)\.bam$"

# generate slopped bed file
echo "slopping bed around patient-specific loci by ${SLOP_BP} bp"
$BEDTOOLS slop -i $BED -g ${HG19_GENOME} -b ${SLOP_BP} | perl -p -i -e 's/ /\t/g' > bed/temp.slopped.bed

#filter lines with >100 bp in the slopped bed
N_BAD_SLOP_LOCI=$(awk '{if ($3-$2 > 100) print $0}' bed/temp.slopped.bed | wc -l)
echo "removing $N_BAD_SLOP_LOCI loci that are off the chromosome"
awk '{if ($3-$2 < 100) print $0}' bed/temp.slopped.bed > bed/temp.slopped.filtered.bed

SLOPPED_BED="bed/temp.slopped.filtered.bed"

# check that length of slopped bed is the same as the original bed
if [ "$(wc -l $BED | awk '{print $1}' - )" != "$(($(wc -l $SLOPPED_BED | awk '{print $1}' -) + $N_BAD_SLOP_LOCI))" ]
  then
  echo "bed length not equal"
  exit 1
fi

while read FILE; do
    if [[ $FILE =~ $regex ]]
    then
        PREFIX="${BASH_REMATCH[1]}"
	else
		echo "${FILE} is not a bam"
		exit 1;
	fi

	# remove paths from prefix for output
	NAME=$(basename $PREFIX)
  echo "NAME = $NAME"

  #code to handle multiple panels (not needed for most projects)
  if [  -f "bed/SLX_barcode.panel.txt" ]
   then
     echo "multiple panels used"

     SLX_BARCODE=$(echo $NAME | grep -o "SLX-....._D..._D..." - )

     PANEL=$(grep $SLX_BARCODE bed/SLX_barcode.panel.txt | awk '{print $2}' -)
     echo "panel used: $PANEL"

     if [ -z "$PANEL" ]
      then
        echo "no panel for this sample"
        exit 1
      fi

     # now grep out that panel from the BED
     echo "cat bed/temp.slopped.bed | grep $PANEL > bed/temp.$PANEL.slopped.bed"
     cat bed/temp.slopped.bed | grep $PANEL > bed/temp.$PANEL.slopped.bed

     SLOPPED_BED="bed/temp.$PANEL.slopped.bed"
     echo "slopped bed (multiple panels) = $SLOPPED_BED"

     NAME=$NAME.$PANEL
   else
    echo "bed file = $BED"
   fi

	echo "file = $FILE"
  echo "file prefix including path = $PREFIX"
  echo "file prefix excluding path = $NAME"
  echo "MQ = $MAPQ"
  echo "BQ = $BASEQ"

  # check that the file has an SLX
  SLX=$(echo $FILE | grep -oP "SLX-[0-9]+")

  if [[ -z $SLX ]] ; then
    echo "ERROR - Files are not correctly named - they need to have SLX-XXXX_DXXX_DXXX"
    #exit
  fi

	#bsub -o pileups.out -e pileups.err $MAKE_PILEUP $PREFIX $BED $OUT $NAME;
  echo "sbatch -o $NAME.pileups.out -e $NAME.pileups.err $MAKE_PILEUP $PREFIX $SLOPPED_BED $OUT $NAME $MAPQ $BASEQ"
  sbatch -o temp/$NAME.pileups.out -e temp/$NAME.pileups.err --job-name="INVAR1" --mem=4G --time="1-1" $MAKE_PILEUP $PREFIX $SLOPPED_BED $OUT $NAME $MAPQ $BASEQ;
done < $INPUT_FILES
