#########################################################
####                  Define output                  ####
#########################################################
# mutect2_joint_calling_output = output+"mutect2_joint_calling/"

#########################################################
####  Define workflow start, stop and error actions  ####
#########################################################
onstart:
    shell("mkdir -p "+mutect2_joint_calling_output+"split/{{vcf,bam,f1r2}}_split")
    shell("mkdir -p "+mutect2_joint_calling_output+"vcf_{{files,filterFlag,filterFlag_PASS,filterFlag_norm_decomposed}}")
    # shell("mkdir -p "+mutect2_joint_calling_output+"stats/")
    shell("mkdir -p "+mutect2_joint_calling_output+"realigned_bam/")


#########################################################
####                  Run All Rules                  ####
#########################################################
# rule all_pairs:
#     input:
#         expand(mutect2_joint_calling_output+"vcf_filterFlag_PASS/{sample}_somatic_mutect2_filterFlag_PASSonly.vcf.gz", sample=SAMPLES)
        # expand(mutect2_joint_calling_output+"vcf_filterFlag_norm_decomposed/{sample}_somatic_mutect2_filterFlag_norm.vcf.gz", sample=SAMPLES)
        # expand(mutect2_joint_calling_output+"vcf_files/{sample}_somatic_mutect2.vcf.gz", sample=SAMPLES),
        # expand(mutect2_joint_calling_output+"{sample}_somatic_mutect2_filterFlag.vcf.gz", sample=SAMPLES),
        # expand(mutect2_joint_calling_output+"{sample}_somatic_mutect2.bam", sample=SAMPLES)
        # expand(mutect2_joint_calling_output+"{sample}_somatic_mutect2.vcf", sample=SAMPLES),
        # expand(mutect2_joint_calling_output+"{sample}_merged_contamination.table", sample=SAMPLES)


##########################################################
####  Call Somatic Variants using Mutect2 on matched  ####
####   Tumor-Normal samples on per chromesome basis   ####
##########################################################
# Assumes PON already calculated in other mutect2 script.
'''
Mutect2 on matched Tumor-Normal samples
'''
rule Mutect2_joint_calling:
    input:
        normal=lambda wildcards: expand("bam/{normal}.connor.recalibrated.bam", normal=config[wildcards.sample]["normal"]),
        tumors=lambda wildcards: expand("bam/{tumor}.connor.recalibrated.bam", tumor=config[wildcards.sample]["tumor"]),
        # pon=pon_location+"pon.vcf.gz",
        intervals=target_regions+"{interval}.interval_list"
    output:
        vcf=mutect2_joint_calling_output+"split/vcf_split/{sample}_somatic_mutect2__{interval}__split.vcf.gz",
        idx=mutect2_joint_calling_output+"split/vcf_split/{sample}_somatic_mutect2__{interval}__split.vcf.gz.tbi",
        bam_subfile=mutect2_joint_calling_output+"split/bam_split/{sample}__{interval}_mutect2.bam",
        vcf_stats=mutect2_joint_calling_output+"split/vcf_split/{sample}_somatic_mutect2__{interval}__split.vcf.gz.stats",
        f1r2=mutect2_joint_calling_output+"split/f1r2_split/{sample}_f1r2__{interval}__split.tar.gz"
    params:
        tumors=lambda wildcards: expand("-I bam/{tumor}.connor.recalibrated.bam", tumor=config[wildcards.sample]["tumor"]),
        normal_name=lambda wildcards: expand("-normal {normal}", normal=config[wildcards.sample]["normal"])
    resources: cpus=3, mem_mb=18000
    shell:
        """
        gatk --java-options -Xmx{resources.mem_mb}g Mutect2 \
        -R {ref} \
        {params.tumors} \
        -I {input.normal} \
        {params.normal_name} \
        -pon {pon_1000g} \
        --germline-resource {gnomad} \
        --af-of-alleles-not-in-resource 0.0000025 \
        --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --native-pair-hmm-threads {resources.cpus} \
        --f1r2-tar-gz {output.f1r2} \
        -L {input.intervals} \
        -bamout {output.bam_subfile} \
        -O {output.vcf}
        """
        # -bamout {output.bam_subfile} \ ## Debugging
        # --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
        # --af-of-alleles-not-in-resource 0.0000025 \
        # --max_alt_alleles_in_normal_count 1000000 \ # NOT IN GATK 4.1 (or 4.2)
        # --max_alt_allele_in_normal_fraction 0.10 \ # NOT IN GATK 4.1 (or 4.2)
        # --tumor-lod-to-emit  \

# rule merge_somatic_vcf:
#     input:
#         vcf_subfile=expand(mutect2_joint_calling_output+"split/vcf_split/{{sample}}_somatic_mutect2__{interval}__split.vcf.gz", interval=INTERVALS)
#     output:
#         vcf=mutect2_joint_calling_output+"vcf_files/{sample}_somatic_mutect2.vcf.gz",
#         idx=mutect2_joint_calling_output+"vcf_files/{sample}_somatic_mutect2.vcf.gz.tbi",
#     params:
#         vcf_subfile=expand("-I "+mutect2_joint_calling_output+"split/vcf_split/{{sample}}_somatic_mutect2__{interval}__split.vcf.gz", interval=INTERVALS)
#     shell:
#         """
#         gatk --java-options {mem} GatherVcfs \
#         {params.vcf_subfile} \
#         -O {output.vcf}

#         tabix -p vcf {output.vcf}
#         """

# ## FOR DEBUGGING
# rule merge_mutect2_bam:
#   input:
#         bam_subfile=expand(mutect2_joint_calling_output+"split/bam_split/{{sample}}__{interval}_mutect2.bam", interval=INTERVALS)
#   output:
#         bam=mutect2_joint_calling_output+"realigned_bam/{sample}_somatic_mutect2.bam"
#   params:
#         bam_subfile=expand("-I "+mutect2_joint_calling_output+"split/bam_split/{{sample}}__{interval}_mutect2.bam", interval=INTERVALS)
#   shell:
#         """
#         gatk --java-options {mem} MergeSamFiles \
#         {params.bam_subfile} \
#         --VALIDATION_STRINGENCY LENIENT \
#         --USE_THREADING true \
#         -O {output.bam} \
#         --CREATE_INDEX true
#         """

# rule merge_somatic_vcf_stats:
#     input:
#         vcf_subfile=expand(mutect2_joint_calling_output+"split/vcf_split/{{sample}}_somatic_mutect2__{interval}__split.vcf.gz.stats", interval=INTERVALS)
#     output:
#         vcf_stats=mutect2_joint_calling_output+"vcf_files/{sample}_somatic_mutect2.vcf.gz.stats",
#     params:
#         vcf_subfile=expand("-stats "+mutect2_joint_calling_output+"split/vcf_split/{{sample}}_somatic_mutect2__{interval}__split.vcf.gz.stats", interval=INTERVALS)
#     shell:
#         """
#         gatk --java-options {mem} MergeMutectStats \
#         {params.vcf_subfile} \
#         -O {output.vcf_stats}
#         """

# #########################################################
# ####          Learn Read Orientation Bias            ####
# ####                                                 ####
# #### Note: Used to fix orientation bias artifacts    ####
# ####       from formalin Formalin-Fixed Paraffin-    ####
# ####       Embedded (FFPE) samples - i.e. not needed ####
# ####       for frozen tissue                         ####
# #########################################################
# rule learnReadOrientationModel:
#     input:
#         f1r2=expand(mutect2_joint_calling_output+"split/f1r2_split/{{sample}}_f1r2__{interval}__split.tar.gz", interval=INTERVALS)
#     output:
#         f1r2_model=mutect2_joint_calling_output+"read_orientation_models/{sample}_read-orientation-model.tar.gz"
#     params:
#         f1r2=expand("-I "+mutect2_joint_calling_output+"split/f1r2_split/{{sample}}_f1r2__{interval}__split.tar.gz", interval=INTERVALS)
#     shell:
#         """
#         gatk LearnReadOrientationModel \
#         {params.f1r2} \
#         -O {output}
#         """
#         # -I {input} \

# #########################################################
# ####           Create Contamination table            ####
# #########################################################
# ### Sæt sammen til en metode, som tager højde for med og uden normal... (Se gatk4 scripts)
# rule GetPileupSummaries:
#     input:
#         bam="bam/{sample}.connor.recalibrated.bam"
#     output:
#         pileup=mutect2_joint_calling_output+"pileups/{sample}_pileup.table"
#     shell:
#         """
#         gatk --java-options {mem} GetPileupSummaries \
#         -I {input.bam} \
#         -V {common_variants} \
#         -L {common_variants} \
#         -O {output}
#         """

# rule CalculateContamination:
#     input:
#         normal=mutect2_joint_calling_output+"pileups/{normal}_pileup.table",
#         tumor=mutect2_joint_calling_output+"pileups/{tumor}_pileup.table"
#     output:
#         contamination=mutect2_joint_calling_output+"contaminations_tables/{tumor}_vs_{normal}_contamination.table"
#     shell:
#         """
#         gatk --java-options {mem} CalculateContamination \
#         -I {input.tumor} \
#         -matched {input.normal} \
#         -O {output}
#         """

# rule JoinContaminationTables:
#     input:
#         cont_tables=lambda wildcards: expand(mutect2_joint_calling_output+"contaminations_tables/{tumor}_vs_{normal}_contamination.table", normal=config[wildcards.sample]["normal"], tumor=config[wildcards.sample]["tumor"])
#     output:
#         joined_table=mutect2_joint_calling_output+"contaminations_tables/{sample}_merged_contamination.table"
#     # params:
#     #     sample_name=expand('{sample}', sample=SAMPLES)
#     shell:
#         """
#         echo -e 'sample\tcontamination\terror' > {output};
#         awk -F'\t' 'FNR == 2' {input} >> {output}
#         """
#         # awk -F'\t' 'FNR == 2' {params.sample_name}*contamination.table >> {output}


# #########################################################
# ####              Filter Mutect2 Calls               ####
# #########################################################
# rule FilterMutectCalls:
#     input:
#         vcf=mutect2_joint_calling_output+"vcf_files/{sample}_somatic_mutect2.vcf.gz",
#         idx=mutect2_joint_calling_output+"vcf_files/{sample}_somatic_mutect2.vcf.gz.tbi",
#         bam=mutect2_joint_calling_output+"realigned_bam/{sample}_somatic_mutect2.bam",
#         stats=mutect2_joint_calling_output+"vcf_files/{sample}_somatic_mutect2.vcf.gz.stats",
#         contamination=mutect2_joint_calling_output+"contaminations_tables/{sample}_merged_contamination.table",
#         read_orientation=mutect2_joint_calling_output+"read_orientation_models/{sample}_read-orientation-model.tar.gz"
#     output:
#         vcf=mutect2_joint_calling_output+"vcf_filterFlag/{sample}_somatic_mutect2_filterFlag.vcf.gz",
#         tbi=mutect2_joint_calling_output+"vcf_filterFlag/{sample}_somatic_mutect2_filterFlag.vcf.gz.tbi",
#         # vcf_pass=mutect2_joint_calling_output+"vcf_filterFlag_PASS/{sample}_somatic_mutect2_filterFlag_PASSonly.vcf.gz",
#         # tbi_pass=mutect2_joint_calling_output+"vcf_filterFlag_PASS/{sample}_somatic_mutect2_filterFlag_PASSonly.vcf.gz.tbi"
#     shell:
#         """
#         gatk --java-options {mem} FilterMutectCalls \
#         -R {ref} \
#         -V {input.vcf} \
#         --contamination-table {input.contamination} \
#         --orientation-bias-artifact-priors {input.read_orientation} \
#         --stats {input.stats} \
#         -L {interval_list} \
#         -O {output.vcf}
#         """
#         # zcat {output.vcf} | egrep '^#|PASS' | bgzip -c > {output.vcf_pass}
#         # tabix -p vcf {output.vcf_pass}

# rule vcf_pass_only:
#     input:
#         vcf=mutect2_joint_calling_output+"vcf_filterFlag/{sample}_somatic_mutect2_filterFlag.vcf.gz",
#     output:
#         vcf_pass=mutect2_joint_calling_output+"vcf_filterFlag_PASS/{sample}_somatic_mutect2_filterFlag_PASSonly.vcf.gz",
#     shell:
#         """
#         bcftools view -i "%FILTER='PASS' | %FILTER='.'" {input.vcf} | bgzip -c > {output.vcf_pass}
#         tabix -p vcf {output.vcf_pass}
#         """

# rule vcf_norm_decompose:
#     input:
#         vcf=mutect2_joint_calling_output+"vcf_filterFlag/{sample}_somatic_mutect2_filterFlag.vcf.gz"
#     output:
#         norm_vcf=mutect2_joint_calling_output+"vcf_filterFlag_norm_decomposed/{sample}_somatic_mutect2_filterFlag_norm.vcf.gz"
#     shell:
#         """
#         zless {input.vcf} \
#         | vt decompose -s - \
#         | vt normalize -r {ref} - \
#         | bgzip -c > {output.norm_vcf}
#         """
#         # vcf index: tabix -p vcf {output.norm_vcf}
#         # Alternative to vt
#         # bcftools norm -m-any {input.vcf} > {output.vcf}

# rule vcf_norm_pass_only:
#     input:
#         norm_vcf=mutect2_joint_calling_output+"vcf_filterFlag_norm_decomposed/{sample}_somatic_mutect2_filterFlag_norm.vcf.gz"
#     output:
#         vcf_pass=mutect2_joint_calling_output+"vcf_filterFlag_norm_decomposed/{sample}_somatic_mutect2_filterFlag_norm_PASSonly.vcf.gz"
#     shell:
#         """
#         bcftools view -i "%FILTER='PASS' | %FILTER='.'" {input.norm_vcf} | bgzip -c > {output.vcf_pass}
#         tabix -p vcf {output.vcf_pass}
#         """