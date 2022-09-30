#!/bin/bash
#===================================================================================
#
# FILE: INVAR3.sh
# USAGE: ./INVAR3.sh
#
# DESCRIPTION: Run Rscript for blacklisting loci and splitting data into
#         target and off target.
#
# REQUIREMENTS: You must check the filtering parameters of the R/INVAR3.R script.
#               The default filtering settings may not be suitable for your data.
#
#===================================================================================

# load config.sh
echo "loading config.sh"
source config/config.sh

mkdir temp
mkdir output_R

sbatch -o output_R/INVAR3.out -e output_R/INVAR3.err --job-name="INVAR3" --time="1-1" --mem=100G --wrap "Rscript $INVAR3_R $SCRIPTS_FOLDER $FINAL_PREFIX $BED $SLX_TABLE"

echo -e "`tput setaf 1``tput bold`you must read output_R/INVAR3.out to see the effects of the filters\n`tput sgr 0`"
