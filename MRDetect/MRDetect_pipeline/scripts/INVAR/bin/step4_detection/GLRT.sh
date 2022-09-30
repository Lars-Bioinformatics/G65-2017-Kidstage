#!/bin/bash
#===================================================================================
#
# FILE: GLRT.sh
# DESCRIPTION: Run Eyal's GLRT R code to calculate likelihood ratios for ctDNA detection
#              This will only run after OUTLIER_SUPPRESSION.sh has first completed
#
# COLLECTING RESULTS:
# `source config/config.sh; DATE=$(date "+%y%m%d"); grep "INVAR_SCORE" GLRT/* > ${FINAL_PREFIX}.INVAR_scores.${DATE}.txt`
#==================================================================================

echo "loading config.sh"
source config/config.sh

echo "beginning of GLRT"

SIZE_COMBINED_PREFIX=$(echo ${SIZE_COMBINED} |sed 's/.rds//g' - | rev | cut -d "/" -f 1 | rev | sed 's/.size_characterisation//g' -)

echo "SMOOTHING = $SMOOTH; sleep for 3sec"

if [ -z "$SMOOTH" ]
  then
    echo "Size smoothing factor not present - please check"
    exit 1
fi

sleep 3s

rm -rf GLRT
mkdir GLRT

rm -rf output_R/downsampled_mixture
mkdir output_R/downsampled_mixture

for FILE in $(ls output_R/os/*os.Rdata); do
  FILE_PREFIX=$(echo $FILE |sed 's/.Rdata//g' - | rev | cut -d "/" -f 1 | rev)

  echo "sbatch -o GLRT/$FILE_PREFIX.GLRT.out -e GLRT$FILE_PREFIX.GLRT.err --wrap Rscript $GLRT_R $FILE $SIZE_COMBINED $DETECTION_FUNCTIONS_R $SMOOTH ${SIZE_COMBINED_PREFIX} ${FINAL_PREFIX}"
  sbatch -o GLRT/$FILE_PREFIX.GLRT.out -e GLRT/$FILE_PREFIX.GLRT.err --job-name="GLRT" --mem=16G --time="1-1" --wrap "Rscript $GLRT_R $FILE $SIZE_COMBINED $DETECTION_FUNCTIONS_R $SMOOTH ${SIZE_COMBINED_PREFIX} ${FINAL_PREFIX}"
done
