#########################################################
####                  Define output                  ####
#########################################################


#########################################################
####  Define workflow start, stop and error actions  ####
#########################################################
onstart:
    shell("mkdir -p "+mutect2_matched_output+"split/{{vcf,bam,f1r2}}_split")
    shell("mkdir -p "+mutect2_matched_output+"vcf_{{files,filterFlag,filterFlag_PASS,filterFlag_norm}}")
    shell("mkdir -p "+mutect2_matched_output+"realigned_bam/")


#########################################################
####                  Run All Rules                  ####
#########################################################
'''
Rule all
'''
# rule all_pairs:
#     input:
#         [expand(mutect2_matched_output+"vcf_filterFlag_PASS/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_PASSonly.vcf.gz",
#             normal=config[fam]["normal"],
#             tumor=config[fam]["tumor"]) for fam in PAIR]


##########################################################
####  Call Somatic Variants using Mutect2 on matched  ####
####   Tumor-Normal samples on per chromesome basis   ####
##########################################################
'''
Mutect2 on matched Tumor-Normal samples
'''
rule Mutect2_matched:
    input:
        normal="bam/{normal}.connor.recalibrated.bam",
        tumor="bam/{tumor}.connor.recalibrated.bam",
        intervals=target_regions+"{interval}.interval_list"
    output:
        vcf=mutect2_matched_output+"split/vcf_split/{tumor}_vs_{normal}_somatic_mutect2__{interval}__split.vcf.gz",
        idx=mutect2_matched_output+"split/vcf_split/{tumor}_vs_{normal}_somatic_mutect2__{interval}__split.vcf.gz.tbi",
        vcf_stats=mutect2_matched_output+"split/vcf_split/{tumor}_vs_{normal}_somatic_mutect2__{interval}__split.vcf.gz.stats",
        bam_subfile=mutect2_matched_output+"split/bam_split/{tumor}_vs_{normal}__{interval}_mutect2.bam",
        f1r2=mutect2_matched_output+"split/f1r2_split/{tumor}_vs_{normal}_f1r2__{interval}__split.tar.gz"
    resources: cpus=3, mem_mb=18000
    shell:
        """
        gatk --java-options -Xmx{resources.mem_mb}g Mutect2 \
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
        -bamout {output.bam_subfile} \
        -L {input.intervals} \
        -O {output.vcf}
        """



# rule computeMatched:
#     input:
#         lambda wildcards: expand(mutect2_matched_output+"vcf_filterFlag_PASS/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_PASSonly.vcf.gz",
#                tumor=config[wildcards.sample]["tumor"],
#                normal=config[wildcards.sample]["normal"])
#     output:
#         "{output}mutect2_{sample}_finished.txt"
#     shell:
#         """
#         echo "SUCCESS" > {output}
#         """



# rule Mutect2_matched:
#     input:
#         normal="bam/{normal}.connor.recalibrated.bam",
#         tumor="bam/{tumor}.connor.recalibrated.bam",
#         intervals=target_regions+"{interval}.interval_list"
#     output:
#         vcf=mutect2_matched_output+"split/vcf_split/{tumor}_vs_{normal}_somatic_mutect2__{interval}__split.vcf.gz",
#         idx=mutect2_matched_output+"split/vcf_split/{tumor}_vs_{normal}_somatic_mutect2__{interval}__split.vcf.gz.tbi",
#         vcf_stats=mutect2_matched_output+"split/vcf_split/{tumor}_vs_{normal}_somatic_mutect2__{interval}__split.vcf.gz.stats",
#         bam_subfile=mutect2_matched_output+"split/bam_split/{tumor}_vs_{normal}__{interval}_mutect2.bam",
#         f1r2=mutect2_matched_output+"split/f1r2_split/{tumor}_vs_{normal}_f1r2__{interval}__split.tar.gz"
#     resources: cpus=3, mem_mb=18000
#     shell:
#         """
#         gatk --java-options -Xmx{resources.mem_mb}g Mutect2 \
#         -R {ref} \
#         -I {input.tumor} \
#         -I {input.normal} \
#         -normal {wildcards.normal} \
#         -pon {pon_1000g} \
#         --germline-resource {gnomad} \
#         --af-of-alleles-not-in-resource 0.0000025 \
#         --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
#         --native-pair-hmm-threads {resources.cpus} \
#         --f1r2-tar-gz {output.f1r2} \
#         -bamout {output.bam_subfile} \
#         -L {input.intervals} \
#         -O {output.vcf}
#         """

# rule merge_somatic_vcf:
#     input:
#         vcf_subfile=expand(mutect2_matched_output+"split/vcf_split/{{tumor}}_vs_{{normal}}_somatic_mutect2__{interval}__split.vcf.gz", interval=INTERVALS)
#     output:
#         vcf=mutect2_matched_output+"vcf_files/{tumor}_vs_{normal}_somatic_mutect2.vcf.gz",
#         idx=mutect2_matched_output+"vcf_files/{tumor}_vs_{normal}_somatic_mutect2.vcf.gz.tbi",
#     params:
#         vcf_subfile=expand("-I "+mutect2_matched_output+"split/vcf_split/{{tumor}}_vs_{{normal}}_somatic_mutect2__{interval}__split.vcf.gz", interval=INTERVALS)
#     shell:
#         """
#         gatk --java-options {mem} GatherVcfs \
#         {params.vcf_subfile} \
#         -O {output.vcf}

#         tabix -p vcf {output.vcf}
#         """

# rule merge_mutect2_bam:
#   input:
#         bam_subfile=expand(mutect2_matched_output+"split/bam_split/{{tumor}}_vs_{{normal}}__{interval}_mutect2.bam", interval=INTERVALS)
#   output:
#         bam=mutect2_matched_output+"realigned_bam/{tumor}_vs_{normal}_somatic_mutect2.bam"
#   params:
#         bam_subfile=expand("-I " + mutect2_matched_output + "split/bam_split/{{tumor}}_vs_{{normal}}__{interval}_mutect2.bam", interval=INTERVALS)
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
#         vcf_subfile=expand(mutect2_matched_output+"split/vcf_split/{{tumor}}_vs_{{normal}}_somatic_mutect2__{interval}__split.vcf.gz.stats", interval=INTERVALS)
#     output:
#         vcf_stats=mutect2_matched_output+"vcf_files/{tumor}_vs_{normal}_somatic_mutect2.vcf.gz.stats",
#     params:
#         vcf_subfile=expand("-stats "+mutect2_matched_output+"split/vcf_split/{{tumor}}_vs_{{normal}}_somatic_mutect2__{interval}__split.vcf.gz.stats", interval=INTERVALS)
#     shell:
#         """
#         gatk --java-options {mem} MergeMutectStats \
#         {params.vcf_subfile} \
#         -O {output.vcf_stats}
#         """

#########################################################
####          Learn Read Orientation Bias            ####
####                                                 ####
#### Note: Used to fix orientation bias artifacts    ####
####       from Formalin-Fixed Paraffin-Embedded     ####
####       (FFPE) samples - i.e. not needed for      ####
####       frozen tissue                             ####
#########################################################
# rule learnReadOrientationModel:
#     input:
#         f1r2=expand(mutect2_matched_output+"split/f1r2_split/{{tumor}}_vs_{{normal}}_f1r2__{interval}__split.tar.gz", interval=INTERVALS)
#     output:
#         f1r2_model=mutect2_matched_output+"read_orientation_models/{tumor}_vs_{normal}_read-orientation-model.tar.gz"
#     params:
#         f1r2=expand("-I "+mutect2_matched_output+"split/f1r2_split/{{tumor}}_vs_{{normal}}_f1r2__{interval}__split.tar.gz", interval=INTERVALS)
#     shell:
#         """
#         gatk LearnReadOrientationModel \
#         {params.f1r2} \
#         -O {output}
#         """

#########################################################
####           Create Contamination table            ####
#########################################################
# rule GetPileupSummaries:
#     input:
#         bam="bam/{sample}.connor.recalibrated.bam"
#     output:
#         pileup=mutect2_matched_output+"pileups/{sample}_pileup.table"
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
#         normal=mutect2_matched_output+"pileups/{normal}_pileup.table",
#         tumor=mutect2_matched_output+"pileups/{tumor}_pileup.table"
#     output:
#         contamination=mutect2_matched_output+"contaminations_tables/{tumor}_vs_{normal}_contamination.table"
#     shell:
#         """
#         gatk --java-options {mem} CalculateContamination \
#         -I {input.tumor} \
#         -matched {input.normal} \
#         -O {output}
#         """

# rule JoinContaminationTables:
#     input:
#         cont_tables=mutect2_matched_output+"contaminations_tables/{tumor}_vs_{normal}_contamination.table"
#     output:
#         joined_table=mutect2_matched_output+"contaminations_tables/{tumor}_vs_{normal}_merged_contamination.table"
#     shell:
#         """
#         echo -e 'sample\tcontamination\terror' > {output};
#         awk -F'\t' 'FNR == 2' {input} >> {output}
#         """


#########################################################
####              Filter Mutect2 Calls               ####
#########################################################
# rule FilterMutectCalls:
#     input:
#         vcf=mutect2_matched_output+"vcf_files/{tumor}_vs_{normal}_somatic_mutect2.vcf.gz",
#         idx=mutect2_matched_output+"vcf_files/{tumor}_vs_{normal}_somatic_mutect2.vcf.gz.tbi",
#         bam=mutect2_matched_output+"realigned_bam/{tumor}_vs_{normal}_somatic_mutect2.bam",
#         stats=mutect2_matched_output+"vcf_files/{tumor}_vs_{normal}_somatic_mutect2.vcf.gz.stats",
#         contamination=mutect2_matched_output+"contaminations_tables/{tumor}_vs_{normal}_merged_contamination.table",
#         read_orientation=mutect2_matched_output+"read_orientation_models/{tumor}_vs_{normal}_read-orientation-model.tar.gz"
#     output:
#         vcf=mutect2_matched_output+"vcf_filterFlag/{tumor}_vs_{normal}_somatic_mutect2_filterFlag.vcf.gz",
#         tbi=mutect2_matched_output+"vcf_filterFlag/{tumor}_vs_{normal}_somatic_mutect2_filterFlag.vcf.gz.tbi"
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

# rule vcf_pass_only:
#     input:
#         vcf=mutect2_matched_output+"vcf_filterFlag/{tumor}_vs_{normal}_somatic_mutect2_filterFlag.vcf.gz",
#     output:
#         vcf_pass=mutect2_matched_output+"vcf_filterFlag_PASS/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_PASSonly.vcf.gz",
#     shell:
#         """
#         bcftools view -i "%FILTER='PASS' | %FILTER='.'" {input.vcf} | bgzip -c > {output.vcf_pass}
#         tabix -p vcf {output.vcf_pass}
#         """

# rule vcf_norm_decompose:
#     input:
#         vcf=mutect2_matched_output+"vcf_filterFlag/{tumor}_vs_{normal}_somatic_mutect2_filterFlag.vcf.gz"
#     output:
#         norm_vcf=mutect2_matched_output+"vcf_filterFlag_norm/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm.vcf.gz"
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