#########################################################
####                  vaf-vaf-plots                  ####
#########################################################

output_readcount=output_raw+"readcount/"
output_vafvafplots = output_dir+"vaf-vaf-plots/"
output_vafvafplots_compare = output_vafvafplots+"ensemble/"
shell("mkdir -p "+output_readcount)
shell("mkdir -p "+output_vafvafplots)
shell("mkdir -p "+output_vafvafplots_compare)

#rule normvars:
#    input:
#        output_variants_notfiltered+"{tumor}_vs_{normal}/{variantcalleroutput}"
#    output:
#        output_variants_notfiltered+"{tumor}_vs_{normal}/norm.{variantcalleroutput}"
#    params:
#        env="bcftools"
#    resources:
#        mem_mb=12000,
#        cpus=threadsparallel
#    threads:
#        threadsparallel
#    shell:
#        """
#        set +eu \
#        && PS1=dummy \
#        && . $(conda info --base)/etc/profile.d/conda.sh \
#        && conda activate {params.env} \
#        && echo $CONDA_PREFIX
#        bcftools view {input} -Oz -o {input}.gz
#        bcftools index -t {input}.gz
#        bcftools norm -m-any {input}.gz | bcftools norm -Ov --check-ref w -f {ref} > {output}
#        """
#
#rule vcf2bed:
#    input: 
#        output_variants_notfiltered+"{tumor}_vs_{normal}/norm.{variantcalleroutput}"
#    output:
#        bed=output_variants_notfiltered+"{tumor}_vs_{normal}/{variantcalleroutput}.bed"
#    resources:
#        mem_mb=12000
#    threads:
#        1
#    shell:
#        """
#        convert2bed-typical \
#        --input=vcf \
#        < {input} \
#        | awk '{{print $1,$2,$3}}' \
#        | sort -k1,1V \
#        | sed --expression='s/ /\t/g' \
#        > {output.bed}
#        """
#
#
#rule mergebed:
#    input: 
#        variantcalleroutput=[expand(expand(output_variants_notfiltered+"{{tumor}}_vs_{{normal}}/{variantcalleroutput}.bed", 
#            variantcalleroutput=VARIANTCALLEROUTPUTS),
#            normal=MATCHED_BIOPSY_NORMAL[1][i],
#            tumor=MATCHED_BIOPSY_NORMAL[0][i]) for i in MATCHED_BIOPSY_NORMAL.index]
#    output:
#        output_variants_notfiltered + "{tumor}_vs_{normal}/all_variants.bed"
#    resources:
#        mem_mb=12000
#    threads:
#        1
#    shell:
#        """
#        bedops -m {input} > {output}
#        """
#
#
#
#rule readcount:
#    input:
#        bam = input_dir + "{tumor}" + bamsuffix,
#        bed = output_variants_notfiltered + "{tumor}_vs_{normal}/all_variants.bed"
#    output:
#        rc = output_variants_notfiltered + "{tumor}_vs_{normal}/{tumor}.rc"
#    resources:
#        mem_mb=12000
#    threads:
#        1
#    shell:
#        """
#        scripts/bam-readcount-master/bin/bam-readcount \
#        -f {ref} \
#        -w 1 \
#        -l {input.bed} \
#        -q 10 \
#        -b 30 \
#        {input.bam} \
#        > {output.rc}
#        """
#
#rule decompose:
#    input:
#        output_variants_notfiltered+"{tumor}_vs_{normal}/norm.{variantcalleroutput}"
#    output:
#        output_variants_notfiltered+"{tumor}_vs_{normal}/decomposed.{variantcalleroutput}"
#    resources:
#        mem_mb=12000
#    threads:
#        1
#    shell:
#        """
#        conda run --name vt \
#        vt decompose \
#        -s {input} \
#        -o {output}
#        """
#
#rule annotate:
#    input:
#        vcf = output_variants_notfiltered+"{tumor}_vs_{normal}/decomposed.{variantcalleroutput}",
#        rc = output_variants_notfiltered + "{tumor}_vs_{normal}/{tumor}.rc"
#    output:
#        output_variants_notfiltered+"{tumor}_vs_{normal}/annotated.{variantcalleroutput}"
#    resources:
#        mem_mb=12000
#    threads:
#        1
#    shell:
#        """
#        vcf-readcount-annotator \
#        {input.vcf} \
#        {input.rc} \
#        DNA \
#        -t all \
#        -s TUMOR \
#        -o {output}
#        """
#
#rule vcf2tsv:
#    input:
#        vcf=output_readcount+"{tumor_all}_vs_{normal_all}_annotated.vcf"
#    output:
#        tsv=output_readcount+"{tumor_all}_vs_{normal_all}_annotated.tsv"
#    resources:
#        mem_mb=12000,
#        cpus=1
#    threads:
#        1
#    shell:
#        """
#        conda run --name vcf \
#        vcf2tsv \
#        -n NA \
#        {input.vcf} \
#        > {output.tsv}
#        """
#
#rule vaf_vaf_plots:
#    input:
#        tsv=[expand(output_readcount+"{tumor_all_pr_normal}_vs_{normal_all_pr_normal}_annotated.tsv",
#            normal_all_pr_normal=normal,
#            prefix_all_pr_normal=MATCHED_TUMOR_NORMAL[2][MATCHED_TUMOR_NORMAL[1] == normal].unique(),
#            tumor_all_pr_normal=MATCHED_TUMOR_NORMAL[0][MATCHED_TUMOR_NORMAL[1] == normal]) for normal in list(dict.fromkeys(MATCHED_TUMOR_NORMAL[1]))]
#    output:
#        pdf=output_vafvafplots+"{prefix_all_pr_normal}_vaf-vaf.pdf"
#    resources:
#        mem_mb=12000,
#        cpus=threadsparallel
#    threads:
#        threadsparallel
#    shell:
#        """
#        Rscript scripts/.internal/vaf-vaf-plot.R \
#        {output.pdf} \
#        {input.tsv}
#        """

rule vcf_cp_ensembles:
    input:
        output_dir+"variants_final/{ensemble}/{tumor}_vs_{normal}/Consensus.sSNV.vcf"
    output:
        output_dir+"variants_final/{ensemble}/{tumor}_vs_{normal}/{ensemble}.vcf"
    resources:
        mem_mb=12000
    threads:
        1
    shell:
        """
        cp {input} {output}
        """

rule vaf_vaf_plots_ensembles:
    input:
        expand(output_dir+"variants_final/{ensemble}/{{tumor}}_vs_{{normal}}/{ensemble}.vcf", ensemble=ENSEMBLE)
    output:
        outdir = directory(output_vafvafplots_compare+"{tumor}_vs_{normal}/"),
        isdone = output_vafvafplots_compare+"{tumor}_vs_{normal}/complete"
    resources:
        mem_mb=12000
    threads:
        threadsparallel
    shell:
        """
        mkdir -p {output.outdir}
        conda run --name r Rscript code/vaf-vaf-plot-compare.R \
        {output.outdir} \
        {input}
        echo "TRUE" > {output.isdone}
        """