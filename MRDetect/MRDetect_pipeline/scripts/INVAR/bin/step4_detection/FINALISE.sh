#!/bin/bash
#===================================================================================
#
# FILE: FINALISE.sh
# USAGE: Run after you have completed all steps of the pipeline to collect files to be analysed locally
#
#
#===================================================================================

echo "loading config.sh"
source config/config.sh

mkdir output_final

echo "copying files into final folder"
# copying error rates file
cp output_R/*.error_rates.Rdata output_final/
# copying on target Rdataframe
cp output_R/*combined.os*.rds output_final/
# copying polsihed rdataframe
cp output_R/*.on_target.Rdata output_final/
# copying size data
cp output_R/*.size_characterisation.rds output_final/
# copying INVAR scores
cp *INVAR_scores*.txt output_final/

echo "now download these files to your local computer"


echo -e "`tput setaf 1``tput bold`moving the size annotation file into output_R folder`tput sgr 0`"
