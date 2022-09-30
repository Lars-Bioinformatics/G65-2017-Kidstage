#########################################################
####                     MeanVAF                     ####
#########################################################

shell("mkdir -p "+output_dir+"meanvaf/")

ALL_PLASMA_MEANVAF=[glob.glob("/mnt/d/master/data/plasmanotfiltered/*.bam")]

rule MeanVAF_vcf2bed:
    input:
        glob.glob(output_dir+"variants_final/*/*/*.sSNV.vcf")
    output:
        output_dir+"meanvaf/all.bed"
    params:
        outdir1=output_dir+"meanvaf/",
        outdir2=output_dir+"meanvaf/tmp_bed/"
    resources:
        mem=6000
    threads:
        1
    shell:
        """
        #export PATH=$PATH:/mnt/c/data/workflows/mvcall/scripts/bedops-master/applications/bed/bedops/bin/:/mnt/c/data/workflows/mvcall/scripts/bedops-master/applications/bed/conversion/bin/:/mnt/c/data/workflows/mvcall/scripts/bedops-master/applications/bed/sort-bed/bin/
        #source ~/.bashrc

        mkdir -p {params.outdir1}
        mkdir -p {params.outdir2}

        i=1
        for file in {input}
        do
            cat $file |
            convert2bed \
            -i vcf \
            -o bed |
            awk '{{ print $1,$2,$3 }}' \
            > {params.outdir2}$i.bed
            i=$(expr $i + 1)
        done

        bedops \
        -m $(ls {params.outdir2}*.bed) |
        sort -k1,1V \
        > {output}

        rm {params.outdir2}*
        rmdir {params.outdir2}
        """

rule MeanVAF_readcount:
    input:
        plasma=ALL_PLASMA_MEANVAF,
        bed=output_dir+"meanvaf/all.bed"
    output:
        end=output_dir+"meanvaf/MeanVAF_readcount_complete"
    params:
        env="readcount",
        outdir=output_dir+"meanvaf/"
    resources:
        mem=6000
    threads:
        1
    shell:
        """
        set +eu \
        && PS1=dummy \
        && . $(conda info --base)/etc/profile.d/conda.sh \
        && conda activate {params.env} \
        && echo $CONDA_PREFIX;

        for FILE in {input.plasma}
        do
            bam-readcount \
            -l {input.bed} \
            -w 1 \
            -f {ref} \
            $FILE \
            > {params.outdir}${{FILE##*/}}.rc
            bam-readcount \
            -l {input.bed} \
            -q 10 \
            -b 30 \
            -w 1 \
            -f {ref} \
            $FILE \
            > {params.outdir}${{FILE##*/}}.qual.rc
        done

        echo "SUCCESS" > {output.end}
        """

rule MeanVAF_annotator:
    input:
        rc=[glob.glob(output_dir+"meanvaf/*.rc")],
        vcf=output_dir+"variants_final/{ensemble}/{tumor}_vs_{normal}/Consensus.sSNV.vcf",
        bed=output_dir+"meanvaf/all.bed",
        end=output_dir+"meanvaf/MeanVAF_readcount_complete"
    output:
        end=output_dir+"meanvaf/{ensemble}/{tumor}_vs_{normal}/MeanVAF_annotator_complete"
    params:
        outdir1=output_dir+"meanvaf/{ensemble}/",
        outdir2=output_dir+"meanvaf/{ensemble}/{tumor}_vs_{normal}/"
    resources:
        mem=6000
    threads:
        1
    shell:
        """
        mkdir -p {params.outdir1}
        mkdir -p {params.outdir2}

        for FILE in {input.rc}
        do
            vcf-readcount-annotator \
            -s TUMOR \
            -o {params.outdir2}${{FILE##*/}}.vcf \
            -t snv \
            {input.vcf} \
            $FILE \
            DNA
        done

        echo "SUCCESS" > {output.end}
        """