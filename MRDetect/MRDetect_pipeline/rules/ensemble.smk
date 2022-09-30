#########################################################
####                     Ensemble                    ####
#########################################################

checkpoint normvars:
    input:
        vcf=output_variants_notfiltered+"{normvars_tumor}_vs_{normvars_normal}/{normvars_variantcalleroutput}.gz"
    output:
        vcf=output_variants_notfiltered+"{normvars_tumor}_vs_{normvars_normal}/norm.{normvars_variantcalleroutput}.gz"
    params:
        env="bcftools"
    resources:
        mem_mb=12000,
        cpus=threadsparallel
    threads:
        threadsparallel
    shell:
        """
        conda run --name bcftools \
        bcftools norm -m-any {input.vcf} | bcftools norm -Oz --check-ref w -f {ref} > {output.vcf}
        conda run --name bcftools \
        bcftools index -t {output.vcf}
        """

def ensemble_input(wildcards):
    return sorted(glob.glob(output_variants_notfiltered + wildcards.ensemble_prefix + "*_vs_" + wildcards.ensemble_normal + "/norm." + wildcards.ensemble_variantcalleroutput + ".gz"))

rule ensemble:
    input:
        vcf=ensemble_input
    output:
        vcf=output_variants_merged+"{ensemble_prefix}_vs_{ensemble_normal}/{ensemble_variantcalleroutput}.gz"
    params:
        outdir=output_variants_merged+"{ensemble_prefix}_vs_{ensemble_normal}/"
    resources:
        mem_mb=12000,
        cpus=threadsparallel
    threads:
        threadsparallel
    shell:
        """
        mkdir -p {params.outdir}
        conda run --name bcftools bcftools isec \
        -n +1 \
        -p {output.vcf}_dir \
        -O z \
        {input.vcf}
        conda run --name bcftools bcftools merge \
        --merge all \
        --force-samples \
        -O z \
        {output.vcf}_dir/*.vcf.gz \
        > {output.vcf}
        """
