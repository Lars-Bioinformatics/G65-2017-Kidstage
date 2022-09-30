#!/bin/bash

##########################################
# Script takes as input case-control pairs
# file to generate detection metrics and
# zscores.

# Minita Shah (mshah@nygenome.org)
# New York Genome Center
##########################################

if [ $# -ne 7 ]; then
    echo ""
    echo "Usage: $0 INP_DIR OUT_DIR TN_PAIRS MEDIAN_THRESH PATIENT_NAME INTERVAL_NAME WIN_SIZE" >&2
    echo ""
    exit 1
fi

inp_dir=$1
out_dir=$2
tn=$3
median_thresh=$4
patient_name=$5
interval_name=$6
win_size=$7

script=`realpath $0`
script_dir=`dirname $script`

summary_out_dir=$out_dir/Summary/robustZscore_CNV_medthresh${median_thresh}/${patient_name}_${interval_name}
if [ ! -d $summary_out_dir ]; then
    mkdir -p $summary_out_dir
fi
temp_dir=$out_dir/Summary/robustZscore_CNV_medthresh${median_thresh}/${patient_name}_${interval_name}/temp
if [ ! -d $temp_dir ]; then
    mkdir -p $temp_dir
fi

numLines=`wc -l $tn | cut -f 1 -d ' '`

for (( f=1; f<=$numLines; f++ )) ;
do
    tumor=`awk -v line=$f 'NR==line && !($1~/^#/) {print $1}' $tn`;
    normal=`awk -v line=$f 'NR==line && !($1~/^#/) {print $2}' $tn`;
    tum_cvg=`awk -v line=$f 'NR==line && !($1~/^#/) {print $3}' $tn`;
    norm_cvg=`awk -v line=$f 'NR==line && !($1~/^#/) {print $4}' $tn`;

	if [ ! -z $tumor ]; then
        # Generating detection metrics per sample--pon pair
        echo $script_dir/parsemedian_robustZscore.py -i $inp_dir -t $tumor -n $normal -in $interval_name -o $summary_out_dir -w $win_size -mt $median_thresh -tc $tum_cvg -nc $norm_cvg
        python $script_dir/parsemedian_robustZscore.py -i $inp_dir -t $tumor -n $normal -in $interval_name -o $summary_out_dir -w $win_size -mt $median_thresh -tc $tum_cvg -nc $norm_cvg
	fi
done

# Converting detection metrics to zscore
echo $script_dir/score_to_detection_zscore.py -s $summary_out_dir -t $tn -o "zscore_summary_metrics.txt"
python $script_dir/score_to_detection_zscore.py -s $summary_out_dir -t $tn -o "zscore_summary_metrics.txt"
