#!/bin/bash
#===================================================================================
#
# FILE: INVAR_SIZE_ANN1.sh
# DESCRIPTION: This shell script generates size information for size_annotation
#       of all fragments, to be used in the INVAR algorithm. Uses a modified
#       script from Dineika Chandrananda.
#
#===================================================================================

rm -rf temp_size

mkdir temp_size
mkdir output_size

echo "loading config.sh"
source config/config.sh

if [[ ! -f to_run.txt ]] ; then
    echo 'File "to_run.txt" does not exist, aborting.'
    echo 'This txt file contains a list of paths to BAM files to be run.'
    exit
fi

regex="^(.+)\.bam$"

awk '{print $1,$3,$4,$5}' $BED | sed 's/ /:/g' | sort | uniq > bed/temp.bed
#BED2=$(cat $BED output_R/$FINAL_PREFIX.off_target.mut_loci.bed)

while read FILE; do
    if [[ $FILE =~ $regex ]]
    then
        PREFIX=$(echo $FILE | rev | cut -d"/" -f 1 | rev)
	else
		echo "${FILE} is not a bam"
		exit 1;
	fi

  echo "python2.7 $SCRIPTS_FOLDER/python/get_fragment_size.py $FILE bed/temp.bed output_size/$PREFIX.inserts_for_annotation.txt"
  sbatch -o temp_size/$PREFIX.get_fragment_size.out \
        -e temp_size/$PREFIX.get_fragment_size.err --job-name="SIZE_ANN1"\
        --time="1-1" --mem=4G --wrap="python2.7 $SCRIPTS_FOLDER/python/get_fragment_size.py $FILE bed/temp.bed output_size/$PREFIX.inserts_for_annotation.csv"

done < $INPUT_FILES
