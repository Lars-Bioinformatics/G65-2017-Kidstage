#########################################################
####                     VarDict                     ####
#########################################################

# Initial operations
vardict_output_dir = output_raw+"vardict/"
vardict_output_tmp = output_raw+"vardict/tmp/"
vardict_output_somatic = output_raw+"vardict/somatic_variants/"
onstart:
    shell("mkdir -p "+vardict_output_dir)
    shell("mkdir -p "+vardict_output_tmp)
    shell("mkdir -p "+vardict_output_somatic)

rule VarDict_call:
    input:
        normal=input_dir+"{normal}"+bamsuffix,
        tumor=input_dir+"{tumor}"+bamsuffix
    output:
        tsv=vardict_output_tmp+"{tumor}_vs_{normal}_somatic_vardict_snvs.tsv"
    params:
        env="vardict"
    resources:
        mem_mb=240000,
        cpus=threadsparallel
    threads:
        threadsparallel
    shell:
        """
        set +eu \
        && PS1=dummy \
        && . $(conda info --base)/etc/profile.d/conda.sh \
        && conda activate {params.env} \
        && echo $CONDA_PREFIX;
        vardict-java \
        -G {ref} \
        -f 0.01 \
        -N {wildcards.tumor} \
        -b "{input.tumor}|{input.normal}" \
        -c 1 \
        -S 2 \
        -E 3 \
        -g 4 \
        -Q 40 \
        -r 5 \
        -th {threads} \
        {medexome} > {output.tsv}
        """

rule VarDict_testsomatic:
    input:
        tsv=vardict_output_tmp+"{tumor}_vs_{normal}_somatic_vardict_snvs.tsv"
    output:
        tsv=vardict_output_tmp+"{tumor}_vs_{normal}_somatic_vardict_snvs_filtered.tsv"
    params:
        env="vardict"
    resources:
        mem_mb=12000,
        cpus=1
    threads:
        1
    shell:
        """
        Rscript scripts/VarDict/testsomatic.R {input.tsv} > {output.tsv}
        """

rule VarDict_var2vcf:
    input:
        tsv=vardict_output_tmp+"{tumor}_vs_{normal}_somatic_vardict_snvs_filtered.tsv"
    output:
        vcf=vardict_output_somatic+"{tumor}_vs_{normal}_somatic_vardict_snvs_filtered.vcf"
    params:
        env="vardict"
    resources:
        mem_mb=12000,
        cpus=1
    threads:
        1
    shell:
        """
        set +eu \
        && PS1=dummy \
        && . $(conda info --base)/etc/profile.d/conda.sh \
        && conda activate {params.env} \
        && echo $CONDA_PREFIX;
        var2vcf_paired.pl \
        -N "{wildcards.tumor}|{wildcards.normal}" \
        -f 0.01 < {input.tsv} > {output.vcf}
        """

checkpoint VarDict_copy:
    input:
        snvs=vardict_output_somatic+"{tumor}_vs_{normal}_somatic_vardict_snvs_filtered.vcf"
    output:
        snvs=output_variants_notfiltered+"{tumor}_vs_{normal}/vardict_bothcombined.vcf"
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