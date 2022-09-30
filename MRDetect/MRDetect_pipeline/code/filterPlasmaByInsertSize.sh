#!/bin/bash
# conda activate sambamba
FILES=$(ls /work/files/plasma/original/*.bam)
n=$(printf "%s\n" $FILES | wc -l)
i=1
for FILE in $FILES
do
    clear
    echo "Filtering input files by insert size (>=90 & <=150) ($i of $n)"
    sambamba view --filter "template_length <= 150 and template_length >= 90" --nthreads 64 --show-progress --format bam --output-filename "/work/files/plasma/min90max150/${FILE##*/}" $FILE
    echo "Filtering input files by insert size (<=150) ($i of $n)"
    sambamba view --filter "template_length <= 150" --nthreads 64 --show-progress --format bam --output-filename "/work/files/plasma/max150/${FILE##*/}" $FILE
    i=$(expr $i + 1)
done