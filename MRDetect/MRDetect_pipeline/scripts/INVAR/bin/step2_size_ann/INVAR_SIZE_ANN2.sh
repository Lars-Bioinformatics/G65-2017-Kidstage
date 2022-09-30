#!/bin/bash
#===================================================================================
#
# FILE: INVAR_SIZE_ANN2.sh
# USAGE: ./INVAR_SIZE_ANN2.sh
#
# DESCRIPTION: Annotate the combined.polished dataframe with size information. Only run once polished df
#               has been made.
#
#===================================================================================

echo "loading config.sh"
source config/config.sh

COMBINED_POLISHED=$(ls output_R/*combined*rds)

rm -rf output_R/combined.polished.size

mkdir output_R/combined.polished.size
mkdir temp

echo "starting size annotation"

for CURR in ${COMBINED_POLISHED}; do
  echo "current file:" $CURR
  SUPPRESSION_SETTING=$(echo $CURR | grep -oP "0.0[0-9]+" -)

  if [ -z "$SUPPRESSION_SETTING" ]
    then
      echo "using default suppression setting"
      SUPPRESSION_SETTING=0.05
  fi

  echo "suppression setting : "$SUPPRESSION_SETTING
  echo "Library prep setting : "$LIBRARY_PREP

  for FILE in $(ls output_size/*inserts_for_annotation*csv); do
    if [ -z ${LIBRARY_PREP+x} ]; then
      echo "Don't recognize library prep will assume Rubicon"
      BC=$(echo $FILE |  grep -oP "D.*_D.*\\." | cut -d "." -f 1)
      SLX=$(echo $FILE | grep -oP "SLX-[0-9]+")
      FILE_PREFIX=${SLX}_${BC}

    elif [ $LIBRARY_PREP == "Rubicon" ]; then
      BC=$(echo $FILE |  grep -oP "D.*_D.*\\." | cut -d "." -f 1)
      SLX=$(echo $FILE | grep -oP "SLX-[0-9]+")
      FILE_PREFIX=${SLX}_${BC}

    elif [ $LIBRARY_PREP == "Agilent" ]; then
      FILE_PREFIX=$(echo $FILE | grep -oP "SLX-[0-9]+_SXTHS[A-Z]0[0-9]")
    else
      echo "Library prep not recognised"
      exit 1

    fi

    echo "FILE = " $FILE
    echo "PREFIX = " $FILE_PREFIX
    echo "COMBINED_POLISHED data = " $CURR
    echo "sbatch -o temp/${FILE_PREFIX}.$SUPPRESSION_SETTING.INVAR_SIZE_MERGING.out -e temp/${FILE_PREFIX}.$SUPPRESSION_SETTING.INVAR_SIZE_MERGING.err --wrap Rscript $MERGE_COMBINED_AND_SIZES $CURR $FILE $FINAL_PREFIX $SCRIPTS_FOLDER $SUPPRESSION_SETTING"

    sbatch -o temp/${FILE_PREFIX}.os_$SUPPRESSION_SETTING.INVAR_SIZE_MERGING.out -e temp/${FILE_PREFIX}.os_$SUPPRESSION_SETTING.INVAR_SIZE_MERGING.err --job-name="SIZE_ANN2" --mem=96G --time="1-1" --wrap "Rscript $MERGE_COMBINED_AND_SIZES $CURR $FILE $FINAL_PREFIX $SCRIPTS_FOLDER $SUPPRESSION_SETTING"
  done

done
