#!/bin/bash
#export PATH=$PATH:/mnt/c/data/workflows/mvcall/scripts/bedops-master/applications/bed/bedops/bin/:/mnt/c/data/workflows/mvcall/scripts/bedops-master/applications/bed/conversion/bin/:/mnt/c/data/workflows/mvcall/scripts/bedops-master/applications/bed/sort-bed/bin/
#source ~/.bashrc

i=1
for file in /mnt/c/data/mvcall-out/pancreas/variants_final/somaticseq_mincallers5/PC1-10-op_tumor_tagseq-medexome_vs_PC1-10-blod_normal_tagseq-medexome/Consensus.sSNV.vcf /mnt/c/data/mvcall-out/pancreas/variants_final/somaticseq_mincallers5/PC1-14-op_tumor_tagseq-medexome_vs_PC1-14-blod_normal_tagseq-medexome/Consensus.sSNV.vcf /mnt/c/data/mvcall-out/pancreas/variants_final/somaticseq_mincallers5/PC1-18-op_tumor_tagseq-medexome_vs_PC1-18-blod_normal_tagseq-medexome/Consensus.sSNV.vcf /mnt/c/data/mvcall-out/pancreas/variants_final/mrdetect/PC1-10-op_tumor_tagseq-medexome_vs_PC1-10-blod_normal_tagseq-medexome/Consensus.sSNV.vcf /mnt/c/data/mvcall-out/pancreas/variants_final/mrdetect/PC1-14-op_tumor_tagseq-medexome_vs_PC1-14-blod_normal_tagseq-medexome/Consensus.sSNV.vcf /mnt/c/data/mvcall-out/pancreas/variants_final/mrdetect/PC1-18-op_tumor_tagseq-medexome_vs_PC1-18-blod_normal_tagseq-medexome/Consensus.sSNV.vcf
do
    cat $file |
    convert2bed \
    -i vcf \
    -o bed |
    awk '{ print $1,$2,$3 }' \
    > /mnt/c/data/mvcall-out/pancreas/variants_final/$i.bed
    i=$(expr $i + 1)
done

bedops \
-m $(ls /mnt/c/data/mvcall-out/pancreas/variants_final/*.bed) |
sort -k1,1V \
> /mnt/c/data/mvcall-out/pancreas/variants_final/all.bed