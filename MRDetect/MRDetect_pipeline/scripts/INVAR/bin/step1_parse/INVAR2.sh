#!/bin/bash
#===================================================================================
#
# FILE: INVAR2.sh
# USAGE: Run once all the final.ann.tsv files have been generated in the output folder
#
# DESCRIPTION: Combine all final.ann.tsv files into one tsv.gz
#
# REQUIREMENTS: None - just make sure INVAR1.sh has fully finished running first.
#
#===================================================================================

# load config.sh
echo "loading config.sh"
source config/config.sh

echo -e "`tput setaf 1``tput bold`this script runs quickly - check output_gz to see if it has finished\n`tput sgr 0`"

mkdir temp

sbatch -o temp/INVAR2.out -e temp/INVAR2.err --job-name="INVAR2" --mem=4G --time="1-1" $CAT_AND_GZIP
