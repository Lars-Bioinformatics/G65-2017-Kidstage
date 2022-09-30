#!/bin/bash
#===================================================================================
#
# FILE: OUTLIER_SUPRESSION.sh
# DESCRIPTION: Outlier suppression step - identify loci with signal that are not consistent with the
#              other patient-specific loci of that sample.
#
#===================================================================================

echo "loading config.sh"
source config/config.sh

if [ -z "$OUTLIER_SUPPRESSION" ]
  then
    echo "OUTLIER_SUPPRESSION.R path not set in config, abort"
    exit 1
fi

# check there were no SLURM memory fails in SIZE_ANN2 step
N_FAILS=$(grep "CANCEL" temp/*MERG*err | wc -l)

if [ "$N_FAILS" -gt "0" ];then
  echo "Some of the files did not complete the size_ann2 step, please check the error files then rerun"
  echo "error files are in temp/"
  #exit 1
fi

rm -rf output_R/os
rm -rf output_R/os_nonzero

mkdir output_R/os
mkdir output_R/os_nonzero

for SLX in $(ls output_R/combined.polished.size/ | grep -oP "SLX-[0-9]+" - | sort |uniq); do
  echo "SLX = $SLX"

  for FILE in $(ls output_R/combined.polished.size/*${SLX}*); do
    FILE_PREFIX=$(echo $FILE | sed 's/.Rdata//g' | cut -d "/" -f 3 - )

    echo "FILE = " $FILE
    echo "PREFIX = " $FILE_PREFIX
    echo "sbatch -o temp/${FILE_PREFIX}.os.out -e temp/${FILE_PREFIX}.os.err --job-name="OS" --mem=4G --time=1-1 --wrap Rscript $OUTLIER_SUPPRESSION $FILE $SCRIPTS_FOLDER"

    sbatch -o temp/${FILE_PREFIX}.os.out -e temp/${FILE_PREFIX}.os.err --job-name="OS" --mem=8G --time="1-1" --wrap "Rscript $OUTLIER_SUPPRESSION $FILE $FILE_PREFIX $SCRIPTS_FOLDER"
  done
  echo "sleep 5s"
  sleep 0.5s

done
