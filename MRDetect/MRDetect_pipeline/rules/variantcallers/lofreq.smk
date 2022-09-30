#########################################################
####                     LoFreq                      ####
#########################################################

# Initial operations
lofreq_output_dir = output_raw+"lofreq/"
onstart:
    shell("mkdir -p "+lofreq_output_dir)

rule LoFreq_somatic:
    input:
        normal=input_dir+"{normal}"+bamsuffix,
        tumor=input_dir+"{tumor}"+bamsuffix,
    output:
        vcf_snvs=lofreq_output_dir+"{tumor}_vs_{normal}/somatic_final_minus-dbsnp.snvs.vcf.gz",
        vcf_indels=lofreq_output_dir+"{tumor}_vs_{normal}/somatic_final_minus-dbsnp.indels.vcf.gz"
    params:
        outdir=lofreq_output_dir+"{tumor}_vs_{normal}/"
    resources:
        mem_mb=12000,
        cpus=1
    threads: 
        1
    shell:
        """
        mkdir -p {params.outdir}
        conda run --name lofreq \
        lofreq somatic \
        -f {ref} \
        -d {dbsnp} \
        -l {medexome} \
        -t {input.tumor} \
        -n {input.normal} \
        -o {params.outdir} \
        --call-indels
        """

# rule LoFreq_somatic:
#     input:
#         normal=input_dir+"{normal}"+bamsuffix,
#         tumor=input_dir+"{tumor}"+bamsuffix,
#         intervals=target_regions+"{interval}.interval_list"
#     output:
#         vcf_snvs=lofreq_output_dir+"split/{tumor}_vs_{normal}/lofreq_final_minus-dbsnp__{interval}__split.snvs.vcf.gz",
#         vcf_indels=lofreq_output_dir+"split/{tumor}_vs_{normal}/lofreq_final_minus-dbsnp__{interval}__split.indels.vcf.gz"
#     params:
#         outdir=lofreq_output_dir+"{tumor}_vs_{normal}/"
#     resources:
#         mem_mb=24000,
#         cpus=1
#     threads: 
#         1
#     shell:
#         """
#         mkdir -p {params.outdir}
#         conda run --name lofreq \
#         lofreq somatic \
#         -f {ref} \
#         -d {dbsnp} \
#         -l {input.intervals} \
#         -t {input.tumor} \
#         -n {input.normal} \
#         -o {params.outdir} \
#         --call-indels
#         """

# rule merge_lofreq_vcf:
#     input:
#         vcf_subfile=expand(lofreq_output_dir+"split/{{tumor}}_vs_{{normal}}/lofreq_final_minus-dbsnp__{interval}__split.{{type}}.vcf.gz", interval=INTERVALS)
#     output:
#         vcf=lofreq_output_dir+"{tumor}_vs_{normal}/lofreq_final_minus-dbsnp.{type}.vcf.gz",
#         tbi=lofreq_output_dir+"{tumor}_vs_{normal}/lofreq_final_minus-dbsnp.{type}.vcf.gz.tbi",
#     params:
#         vcf_subfile=expand("-I "+lofreq_output_dir+"split/{{tumor}}_vs_{{normal}}/lofreq_final_minus-dbsnp__{interval}__split.{{type}}.vcf.gz", interval=INTERVALS)
#     resources:
#         mem_mb=12000,
#         cpus=1
#     shell:
#         """
#         conda run --name gatk4 \
#         gatk --java-options {mem} GatherVcfs \
#         {params.vcf_subfile} \
#         -O {output.vcf}

#         tabix -p vcf {output.vcf}
#         """

checkpoint LoFreq_copy:
    input:
        snvs=lofreq_output_dir+"{tumor}_vs_{normal}/somatic_final_minus-dbsnp.snvs.vcf.gz",
        indels=lofreq_output_dir+"{tumor}_vs_{normal}/somatic_final_minus-dbsnp.indels.vcf.gz",
        # snvs=expand(lofreq_output_dir+"{{tumor}}_vs_{{normal}}/somatic_final_minus-dbsnp.{type}.vcf.gz",type="snvs"),
        # indels=expand(lofreq_output_dir+"{{tumor}}_vs_{{normal}}/somatic_final_minus-dbsnp.{type}.vcf.gz",type="indels")
    output:
        snvs=output_variants_notfiltered+"{tumor}_vs_{normal}/lofreq_snvs.vcf",
        indels=output_variants_notfiltered+"{tumor}_vs_{normal}/lofreq_indels.vcf"
    params:
        folder=output_variants_notfiltered+"{tumor}_vs_{normal}/"
    resources:
        mem_mb=12000,
        cpus=1
    threads:
        1
    shell:
        """
        mkdir -p {params.folder}
        gzip -c -d {input.snvs} > {output.snvs}
        gzip -c -d {input.indels} > {output.indels}
        """
