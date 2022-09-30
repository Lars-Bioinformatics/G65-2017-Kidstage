#!/bin/bash
set +eu \
&& PS1=dummy \
&& . $(conda info --base)/etc/profile.d/conda.sh \
&& conda activate readcount \
&& echo $CONDA_PREFIX;
FILES=$(ls /mnt/d/master/data/plasmafilteredmaxinsertsize/*.bam)
n=$(printf "%s\n" $FILES | wc -l)
i=1
for FILE in $FILES
do
    echo "bam-readcount ($i of $n)"
    bam-readcount \
    -l /mnt/c/data/mvcall-out/pancreas/variants_final/all.bed \
    -w 1 \
    -f /mnt/d/master/resources/Homo_sapiens_assembly38.fasta \
    $FILE \
    > ${FILE%*.}.rc
    i=$(expr $i + 1)
done
i=1
for FILE in $FILES
do
    echo "bam-readcount qual ($i of $n)"
    bam-readcount \
    -l /mnt/c/data/mvcall-out/pancreas/variants_final/all.bed \
    -q 10 \
    -b 30 \
    -w 1 \
    -f /mnt/d/master/resources/Homo_sapiens_assembly38.fasta \
    $FILE \
    > ${FILE%*.}.qual.rc
    i=$(expr $i + 1)
done
echo Success