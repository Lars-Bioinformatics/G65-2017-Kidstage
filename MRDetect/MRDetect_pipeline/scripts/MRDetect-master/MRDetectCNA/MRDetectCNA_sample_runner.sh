#!/bin/bash

###############################################
# Script takes as input CNV bed file, 
# sample bam, reference fasta and bin size to 
# generate binned read counts.

# Minita Shah (mshah@nygenome.org)
# New York Genome Center
##############################################

if [ $# -lt 6 ]; then
    echo ""
    echo "Usage: $0 BAM REF WIN_SIZE CNV_BED OUT_DIR SAMPLE_NAME INTERVAL_NAME" >&2
    echo ""
    exit 1
fi

bam=$1
ref=$2
win_size=$3
cnv_bed=$4
out_dir=$5
sample_name=$6
interval_name=$7

script=`realpath $0`
script_dir=`dirname $script`

sample_dir="$out_dir/$sample_name"
if [ ! -d $sample_dir ]; then
    mkdir -p $sample_dir
fi

temp_dir="$out_dir/$sample_name/temp"
interval_dir="$out_dir/$sample_name/$interval_name/intervals"
segments_out_dir="$out_dir/$sample_name/$interval_name/segments"

if [ ! -d $temp_dir ]; then
    mkdir -p $temp_dir
fi

if [ ! -d $interval_dir ]; then
    mkdir -p $interval_dir
fi

if [ ! -d $segments_out_dir ]; then
    mkdir -p $segments_out_dir
fi

# Convert CNV bed file to intervals
grep -v "^#" $cnv_bed | awk '{print $1,$2,$3,$4}' | while read v1 v2 v3 v4; do cat $script_dir/interval_header.template > $interval_dir/"$v1"_"$v2"_"$v3"."$v4".cnv.interval_list; printf "$v1\t$v2\t$v3\t+\twindow\n" >> $interval_dir/"$v1"_"$v2"_"$v3"."$v4".cnv.interval_list; done

# Run DOC for each CNV interval - running in parallel for each interval
for i in $(ls $interval_dir); do
    interval_n=$(echo $i | awk -F"." '{print $1}')
    interval_type=$(echo $i | awk -F"." '{print $2}')
    if [ ! -f "$segments_out_dir/$interval_n.$interval_type.$s.DepthOfCoverage.sample_summary" ]; then
        #qsub -N $interval_name.$sample_name.doc -l h_rss=32G -o $temp_dir/$interval_n.$sample_name.doc.out.txt -e $temp_dir/$interval_n.$sample_name.doc.err.txt $script_dir/DepthOfCoverage.sh $temp_dir $sample_name $interval_n $bam ${interval_dir}/$i ${interval_type} $segments_out_dir $ref
        sh $script_dir/DepthOfCoverage.sh $temp_dir $sample_name $interval_n $bam ${interval_dir}/$i ${interval_type} $segments_out_dir $ref
    fi
done

# Generate median counts for each window - depends on the output of the previous step
#qsub -hold_jid $interval_name.$sample_name.doc -l h_rss=8G -o $temp_dir/$sample_name.CNV.printmedian.out.txt -e $temp_dir/$sample_name.CNV.printmedian.err.txt $script_dir/printmedian.py -i $segments_out_dir -o $segments_out_dir -s $sample_name -f $cnv_bed -w $win_size -t CNV -z
sh $script_dir/printmedian.py -i $segments_out_dir -o $segments_out_dir -s $sample_name -f $cnv_bed -w $win_size -t CNV -z
#qsub -hold_jid $interval_name.$sample_name.doc -l h_rss=8G -o $temp_dir/$sample_name.NEU.printmedian.out.txt -e $temp_dir/$sample_name.NEU.printmedian.err.txt $script_dir/printmedian.py -i $segments_out_dir -o $segments_out_dir -s $sample_name -f $cnv_bed -w $win_size -t NEU -z
sh $script_dir/printmedian.py -i $segments_out_dir -o $segments_out_dir -s $sample_name -f $cnv_bed -w $win_size -t NEU -z

