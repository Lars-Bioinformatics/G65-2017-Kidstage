#!/bin/bash
#===================================================================================
#
# FILE: INVAR_SIZE_CHARACTERISATION.sh
#
# DESCRIPTION: Generate size_characterisation dataframe from annotated Rdata files.
#              Does not necessarily need to be run if you provide it with a separate
#              size characterisation df.
#
#===================================================================================

echo "loading config.sh"
source config/config.sh

echo "starting size annotation 3 script"

echo "SLX table = " $SLX_TABLE
echo "sbatch -o temp/INVAR_SIZE_CHARACTERISATION.out -e temp/INVAR_SIZE_CHARACTERISATION.err --wrap Rscript $INVAR_SIZE_CHARACTERISATION_R $STUDY_ID $FINAL_PREFIX $SLX_TABLE $SCRIPTS_FOLDER"
sbatch -o temp/INVAR_SIZE_CHARACTERISATION.out -e temp/INVAR_SIZE_CHARACTERISATION.err --job-name="SIZE_ANN3" --time="1-1" --mem=64G --wrap "Rscript $INVAR_SIZE_CHARACTERISATION_R $FINAL_PREFIX $SLX_TABLE $SCRIPTS_FOLDER"
