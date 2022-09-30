#!/bin/bash
#===================================================================================
#
# FILE: ANNOTATE_WITH_OS.sh
# DESCRIPTION: Generate size_characterisation dataframe from annotated Rdata files.
# Does not necessarily need to be run if you provide it with a separate size characterisation df.
#
#===================================================================================

echo "loading config.sh"
source config/config.sh

echo "collecting repolished files"

sbatch -o temp/ANN_OS.out -e temp/ANN_OS.err --job-name="ANN_OS" --time="1-1" --mem=32G --wrap "Rscript $OS_ANNOTATION $SCRIPTS_FOLDER"
