#!/bin/bash

#===================================================================================
#
# FILE: config.sh
#
# DESCRIPTION: This config file contains parameters for running INVAR code. Please fill
#              in the paths below. If you are running this from outside the Rosenfeld lab,
#               please check every path in the config carefully.
#
# REQUIRED SOFTWARE: see invar wiki
#===================================================================================

# Please fill in the variables below [ESSENTIAL]
STUDY_ID="MELR" #
ERROR_SUPPRESION_NAME="f0.9_s2" # unique identifier for CONNOR setting; f = consensus threshold, s = family size
BED="bed/example.bed"  # [bed/example.bed] BED file of patient-specific loci
MAPQ=40 #[40]	min MAPQ threshold
BASEQ=20  #[20] min BASEQ threshold
INPUT_FILES="to_run.txt" # [./to_run.txt] this is a txt file containing the bam IDs to be run. Ensure that all files are in the same folder
LIBRARY_PREP="Rubicon" # choose between Rubicon and Agilent -- needed for drawing out barcode information
SLX_TABLE="combined.SLX_table.180124.csv"
TUMOUR_MUTATIONS_CSV="bed/MELR.patient_specific.baseline.csv"

# Please check over the variables below [non-essential, can use defaults]
CHUNKS=5 #[5] split mpileup output into N files. if you have few files then set this higher
OUT="output" # [output] output folder.
OUTPUT_SIZE_FOLDER="output_size" # [output_size] - output folder for size annotation files
MIN_DP=5 # [5] min DP to consider for mpileup. Set to 1 for sWGS samples.
SLOP_BP=10 #[10] set how many bases either side of the target base to assess for background error rate.
SMOOTH=0.25 # [0.25] smoothing function for size profile - width of smoothing
MAX_DP=100000 # MAX DP for mpileup. MAX DP is also set in config2 when parsing TSV files.

# Do not edit below if you are NRlab ==================================
# set dir of input folder
INPUT_FOLDER=$(head -n 1 ${INPUT_FILES} | rev | cut -d'/' -f2- | rev)

# software paths
MPILEUP="samtools mpileup"
SAMTOOLS_INDEX="samtools index"
BCFTOOLS="bcftools"
PICARD="picard.jar"
BEDTOOLS="bedtools"

# reference paths
FASTAREF="ucsc.hg19.fasta"
HG19_GENOME="hg19.genome"
K1G_DB="ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz" #1000 genomes SNP database
COSMIC_DB="CosmicCodingMuts.vcf.gz" # cosmic database

# Paths within repository to each of the scripts
SCRIPTS_FOLDER="/Users/${USER}/tapas_pipeline" # location of the scripts folder
MAKE_PILEUP="${SCRIPTS_FOLDER}/make_mpileup.sh"
ANNOTATE_SCRIPT="${SCRIPTS_FOLDER}/ann.sh"
ANNOTATE_PYTHON="${SCRIPTS_FOLDER}/python/ann.py"
CAT_AND_GZIP="${SCRIPTS_FOLDER}/cat_and_gzip.sh"
INVAR3_R="${SCRIPTS_FOLDER}/R/INVAR3.R"
ERROR_RATE_R="${SCRIPTS_FOLDER}/R/background_error_rates.R"
MERGE_COMBINED_AND_SIZES="${SCRIPTS_FOLDER}/R/size_annotation.R"
DETECTION_FUNCTIONS_R="${SCRIPTS_FOLDER}/R/detection_functions.R"
GLRT_R="${SCRIPTS_FOLDER}/R/GLRT_iterate_nonptspec.R"
INVAR_SIZE_CHARACTERISATION_R="${SCRIPTS_FOLDER}/R/size_annotation3.R"
OUTLIER_SUPPRESSION="${SCRIPTS_FOLDER}/R/outlier_suppression.R"

# utils - not core parts of the pipeline
OS_ANNOTATION="${SCRIPTS_FOLDER}/R/annotate_with_outlier_suppression.R"
GLRT_DOWNSAMPLED_R="${SCRIPTS_FOLDER}/R/downsample/GLRT.downsample_molecules.R"
GLRT_DOWNSAMPLED_PROPORTIONATE_R="${SCRIPTS_FOLDER}/R/downsample/GLRT.downsample_proportionate.R"
GLRT_DOWNSAMPLED_LOCI_R="${SCRIPTS_FOLDER}/R/downsample/GLRT_downsample_loci.R"

# final file prefix
FINAL_PREFIX="${STUDY_ID}.${ERROR_SUPPRESION_NAME}.BQ_${BASEQ}.MQ_${MAPQ}"
SIZE_COMBINED="output_R/${FINAL_PREFIX}.size_characterisation.rds"

# echo these variables for the user to check
echo "TEXT FILE OF INPUT FILES = " $INPUT_FILES
echo "BED FILE = " $BED
echo "OUTPUT FOLDER = " $OUT
echo "MAPQ = " $MAPQ
echo "BASEQ = " $BASEQ
echo "MIN_DP = " $MIN_DP
echo "ERROR_SUPPRESION_NAME = " $ERROR_SUPPRESION_NAME

# variables useful for 'for loops'
SLEEP_EVERY_N_LINES=${1:-100}
SLEEP_DURATION=${2:-2}

[[ -z "$1" ]] && { echo "Parameter 1 (SLEEP_EVERY_N_LINES) is empty, default to 100" ; }
[[ -z "$2" ]] && { echo "Parameter 1 (SLEEP_DURATION) is empty, default to 2" ; }

# Checking inputs ==================================
# Code below runs checks on input files before running
# check with the user if they wish to run the example data
FIRST_FILE=$(head -n 1 $INPUT_FILES)
echo "first file = " $FIRST_FILE

if test "$FIRST_FILE" = "input/SLX-00000_D000_D000.mrkdup.connor.f0.9_s2.bam"; then
    read -n1 -r -p "Are you sure you wish to run the example files? y/n"  key

    if test "$key" = 'y' ; then
        printf "\nstarting script"
    else
        printf "Please check your input files"
        exit 1
    fi

else
  echo "Running now ..."

fi

# Check that the user has provided the BED for just 1bp, otherwise stop
BED_START_STOP_DIFF=$(tail -n +2 $BED | head -n 1 | awk '{print $3-$2}' -)

if test $BED_START_STOP_DIFF -ne 1; then
        printf "Your BED file is not of 1bp positions, aborting. Check your bed file\n"
        printf "NEEDS TO BE: CHR  POS-1 POS REF ALT"
        printf "No more than 1 line of header! Separated\tby\ttabs!"
        exit 1
fi
