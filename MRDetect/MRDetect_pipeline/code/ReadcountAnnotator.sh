#!/bin/bash
#METHODS=(mrdetect somaticseq_mincallers5)
#SAMPLES=(PC1-10-op_tumor_tagseq-medexome_vs_PC1-10-blod_normal_tagseq-medexome PC1-14-op_tumor_tagseq-medexome_vs_PC1-14-blod_normal_tagseq-medexome PC1-18-op_tumor_tagseq-medexome_vs_PC1-18-blod_normal_tagseq-medexome)
FILES=$(ls /mnt/d/master/data/plasmafilteredmaxinsertsize/*.rc)
n=$(printf "%s\n" $FILES | wc -l)
i=1
for FILE in $FILES
do
    echo "vcf-readcount-annotator ($i of $n)"
    vcf-readcount-annotator \
    -s TUMOR \
    -o /mnt/c/data/mvcall-out/pancreas/variants_final/mrdetect/PC1-10-op_tumor_tagseq-medexome_vs_PC1-10-blod_normal_tagseq-medexome/maxinsertsize150/${FILE##*/}.vcf \
    -t snv \
    /mnt/c/data/mvcall-out/pancreas/variants_final/mrdetect/PC1-10-op_tumor_tagseq-medexome_vs_PC1-10-blod_normal_tagseq-medexome/Consensus.sSNV.vcf \
    $FILE \
    DNA
    i=$(expr $i + 1)
done
