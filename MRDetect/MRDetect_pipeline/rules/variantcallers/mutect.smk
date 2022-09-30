#########################################################
####                     MuTect2                     ####
#########################################################

# Initial operations
mutect2_output_somatic = output_raw+"mutect2/"
onstart:
    shell("mkdir -p "+mutect2_output_somatic)

# Call Somatic Variants using Mutect2 on matched Tumor-Normal samples on per chromesome basis
rule Mutect_matched:
    input:
        normal=input_dir+"{normal}"+bamsuffix,
        tumor=input_dir+"{tumor}"+bamsuffix,
        intervals=target_regions+"{interval}.interval_list"
    output:
        vcf=mutect2_output_somatic+"split/vcf_split/{tumor}_vs_{normal}_somatic_mutect2__{interval}__split.vcf.gz",
        idx=mutect2_output_somatic+"split/vcf_split/{tumor}_vs_{normal}_somatic_mutect2__{interval}__split.vcf.gz.tbi",
        vcf_stats=mutect2_output_somatic+"split/vcf_split/{tumor}_vs_{normal}_somatic_mutect2__{interval}__split.vcf.gz.stats",
        f1r2=mutect2_output_somatic+"split_f1r2/{tumor}_vs_{normal}_f1r2__{interval}__split.tar.gz"
        # bam_subfile=mutect2_output_somatic+"bam_split/{tumor}_vs_{normal}__{interval}_mutect2.bam",
    resources: cpus=4, mem_mb=24000
    shell:
        """
        conda run --name gatk4 \
        gatk --java-options {mem} Mutect2 \
        -R {ref} \
        -I {input.tumor} \
        -I {input.normal} \
        -normal {wildcards.normal} \
        -pon {pon_1000g} \
        --germline-resource {gnomad} \
        --af-of-alleles-not-in-resource 0.0000025 \
        --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --native-pair-hmm-threads {resources.cpus} \
        --f1r2-tar-gz {output.f1r2} \
        -L {input.intervals} \
        -O {output.vcf}
        """
        # -bamout {output.bam_subfile} \ ## Debugging
        # --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
        # --af-of-alleles-not-in-resource 0.0000025 \
        # --max_alt_alleles_in_normal_count 1000000 \ # NOT IN GATK 4.1 (or 4.2)
        # --max_alt_allele_in_normal_fraction 0.10 \ # NOT IN GATK 4.1 (or 4.2)
        # --tumor-lod-to-emit  \

rule merge_somatic_vcf:
    input:
        vcf_subfile=expand(mutect2_output_somatic+"split/vcf_split/{{sample}}_somatic_mutect2__{interval}__split.vcf.gz", interval=INTERVALS)
    output:
        vcf=mutect2_output_somatic+"vcf_files/{sample}_somatic_mutect2.vcf.gz",
        idx=mutect2_output_somatic+"vcf_files/{sample}_somatic_mutect2.vcf.gz.tbi",
    params:
        vcf_subfile=expand("-I "+mutect2_output_somatic+"split/vcf_split/{{sample}}_somatic_mutect2__{interval}__split.vcf.gz", interval=INTERVALS)
    resources:
        mem_mb=6000,
        cpus=1
    shell:
        """
        conda run --name gatk4 \
        gatk --java-options {mem} GatherVcfs \
        {params.vcf_subfile} \
        -O {output.vcf}

        tabix -p vcf {output.vcf}
        """

rule merge_somatic_vcf_stats:
    input:
        vcf_subfile=expand(mutect2_output_somatic+"split/vcf_split/{{sample}}_somatic_mutect2__{interval}__split.vcf.gz.stats", interval=INTERVALS)
    output:
        vcf_stats=mutect2_output_somatic+"vcf_files/{sample}_somatic_mutect2.vcf.gz.stats",
    params:
        vcf_subfile=expand("-stats "+mutect2_output_somatic+"split/vcf_split/{{sample}}_somatic_mutect2__{interval}__split.vcf.gz.stats", interval=INTERVALS)
    resources:
        mem_mb=6000,
        cpus=1
    shell:
        """
        conda run --name gatk4 \
        gatk --java-options {mem} MergeMutectStats \
        {params.vcf_subfile} \
        -O {output.vcf_stats}
        """

# Learn Read Orientation Bias
rule learnReadOrientationModel:
    input:
        f1r2=expand(mutect2_output_somatic+"split_f1r2/{{tumor}}_vs_{{normal}}_f1r2__{interval}__split.tar.gz", interval=INTERVALS)
    output:
        f1r2_model=mutect2_output_somatic+"read_orientation/{tumor}_vs_{normal}_read-orientation-model.tar.gz"
    params:
        f1r2=expand("-I "+mutect2_output_somatic+"split_f1r2/{{tumor}}_vs_{{normal}}_f1r2__{interval}__split.tar.gz", interval=INTERVALS)
    resources:
        mem_mb=6000,
        cpus=1
    shell:
        """
        conda run --name gatk4 \
        gatk LearnReadOrientationModel \
        {params.f1r2} \
        -O {output}
        """


# Create Contamination table
rule Mutect_GetPileupSummaries_normal:
    input:
        bam=input_dir+"{normal}"+bamsuffix
    output:
        pileup=mutect2_output_somatic+"contamination/{normal}_normal_pileup.table"
    resources:
        mem_mb=60000,
        cpus=1
    threads:
        1
    shell:
        """
        conda run --name gatk4 \
        gatk --java-options {mem} GetPileupSummaries \
        -I {input.bam} \
        -V {common_variants} \
        -L {common_variants} \
        -O {output}
        """

rule Mutect_GetPileupSummaries_tumor:
    input:
        bam=input_dir+"{tumor}"+bamsuffix
    output:
        pileup=mutect2_output_somatic+"contamination/{tumor}_tumor_pileup.table"
    resources:
        mem_mb=150000,
        cpus=1
    threads:
        1
    shell:
        """
        conda run --name gatk4 \
        gatk --java-options {mem} GetPileupSummaries \
        -I {input.bam} \
        -V {common_variants} \
        -L {common_variants} \
        -O {output}
        """

rule Mutect_CalculateContamination:
    input:
        normal=mutect2_output_somatic+"contamination/{normal}_normal_pileup.table",
        tumor=mutect2_output_somatic+"contamination/{tumor}_tumor_pileup.table"
    output:
        contamination=mutect2_output_somatic+"contamination/{tumor}_vs_{normal}_contamination.table"
    resources:
        mem_mb=12000,
        cpus=1
    threads:
        1
    shell:
        """
        conda run --name gatk4 \
        gatk --java-options {mem} CalculateContamination \
        -I {input.tumor} \
        -matched {input.normal} \
        -O {output}
        """


# Filter Mutect2 Calls
rule FilterMutectCalls:
    input:
        vcf=mutect2_output_somatic+"vcf_files/{tumor}_vs_{normal}_somatic_mutect2.vcf.gz",
        idx=mutect2_output_somatic+"vcf_files/{tumor}_vs_{normal}_somatic_mutect2.vcf.gz.tbi",
        stats=mutect2_output_somatic+"vcf_files/{tumor}_vs_{normal}_somatic_mutect2.vcf.gz.stats",
        contamination=mutect2_output_somatic+"contamination/{tumor}_vs_{normal}_contamination.table",
        read_orientation=mutect2_output_somatic+"read_orientation/{tumor}_vs_{normal}_read-orientation-model.tar.gz"
    output:
        vcf=mutect2_output_somatic+"{tumor}_vs_{normal}_somatic_mutect2_filtered.vcf.gz",
        idx=mutect2_output_somatic+"{tumor}_vs_{normal}_somatic_mutect2_filtered.vcf.gz.tbi"
    resources:
        mem_mb=12000,
        cpus=1
    shell:
        """
        conda run --name gatk4 \
        gatk --java-options {mem} FilterMutectCalls \
        -R {ref} \
        -V {input.vcf} \
        --contamination-table {input.contamination} \
        --orientation-bias-artifact-priors {input.read_orientation} \
        --stats {input.stats} \
        -L {medexome} \
        -O {output.vcf}
        """

checkpoint Mutect_copy:
    input:
        bothcombined=mutect2_output_somatic+"{tumor}_vs_{normal}_somatic_mutect2_filtered.vcf.gz"
    output:
        bothcombined=output_variants_notfiltered+"{tumor}_vs_{normal}/mutect_bothcombined.vcf"
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
        gzip -c -d {input.bothcombined} > {output.bothcombined}
        # cp {input.bothcombined} {output.bothcombined}
        """