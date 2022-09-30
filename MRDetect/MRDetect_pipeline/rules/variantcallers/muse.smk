#########################################################
####                      MuSE                       ####
#########################################################

# Initial operations
muse_output_dir = output_raw+"muse/"
muse_output_tmp = output_raw+"muse/tmp/"
muse_output_somatic = output_raw+"muse/somatic_variants/"
onstart:
    shell("mkdir -p "+muse_output_dir)
    shell("mkdir -p "+muse_output_tmp)
    shell("mkdir -p "+muse_output_somatic)

rule MuSE_call:
    input:
        normal=input_dir+"{normal}"+bamsuffix,
        tumor=input_dir+"{tumor}"+bamsuffix
    output:
        txt=muse_output_tmp+"{tumor}_vs_{normal}.MuSE.txt"
    params:
        prefix=muse_output_tmp+"{tumor}_vs_{normal}"
    resources:
        mem_mb=240000,
        cpus=1
    threads:
        1
    shell:
        """
        scripts/MuSE/MuSEv1.0rc_submission_c039ffa call \
        -f {ref} \
        -O {params.prefix} \
        {input.tumor} \
        {input.normal}
        """

rule MuSE_sump:
    input:
        f=muse_output_tmp+"{tumor}_vs_{normal}.MuSE.txt"
    output:
        vcf=muse_output_somatic+"{tumor}_vs_{normal}_somatic_muse_snvs.vcf"
    resources:
        mem_mb=24000,
        cpus=1
    threads:
        1
    shell:
        """
        scripts/MuSE/MuSEv1.0rc_submission_c039ffa sump \
        -I {input.f} \
        -E \
        -O {output.vcf} \
        -D {dbsnp}
        """

checkpoint MuSE_copy:
    input:
        snvs=muse_output_somatic+"{tumor}_vs_{normal}_somatic_muse_snvs.vcf"
    output:
        snvs=output_variants_notfiltered+"{tumor}_vs_{normal}/muse_snvs.vcf"
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
        cp {input.snvs} {output.snvs}
        """