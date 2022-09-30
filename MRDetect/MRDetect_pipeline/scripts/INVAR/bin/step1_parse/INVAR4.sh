#!/bin/bash
#===================================================================================
#
# FILE: INVAR4.sh
# USAGE: ./INVAR4.sh
#
# DESCRIPTION: Run Rscript to filter and annotate on_target loci. config2.R needs to be updated
#
#===================================================================================

# load config.sh
echo "loading config.sh"
source config/config.sh

echo "running INVAR4"
echo "sbatch -o output_R/INVAR4.out -e output_R/INVAR4.err --job-name=INVAR4 --wrap Rscript $INVAR4_R $TUMOUR_MUTATIONS_CSV"

sbatch -o output_R/INVAR4.out -e output_R/INVAR4.err --job-name="INVAR4" --time="1-1" --mem=16G --wrap "Rscript ${SCRIPTS_FOLDER}/R/INVAR4.R $TUMOUR_MUTATIONS_CSV $SLX_TABLE"
