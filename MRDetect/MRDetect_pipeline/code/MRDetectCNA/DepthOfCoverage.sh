#!/bin/bash

########################################
# Minita Shah (mshah@nygenome.org)
# New York Genome Center
########################################

if [ $# -ne 8 ]; then
    echo ""
    echo "Usage: $0 TEMP_DIR SAMPLE_NAME SEG_NAME BAM_FILE INTERVAL_LIST INTERVAL_TYPE OUT_DIR REF" >&2
    echo ""
    exit 1
fi

temp_dir=$1
sample_name=$2
seg_name=$3
bam_file=$4
interval_list=$5
interval_type=$6
out_dir=$7
ref=$8

scratch_hash="$temp_dir"

java \
-Djava.io.tmpdir=$scratch_hash \
-Xmx24576m \
-jar /mnt/c/data/workflows/mvcall/scripts/MRDetect-master/MRDetectCNA/gakt3.4/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-dels \
-rf BadCigar \
-ct 100 \
-L $interval_list \
--interval_merging OVERLAPPING_ONLY \
-R $ref \
-I $bam_file \
-o $out_dir/$seg_name.$interval_type.$sample_name.DepthOfCoverage
