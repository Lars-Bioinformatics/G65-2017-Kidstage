__title__ = "Pipeline for Somatic JOINT Variant Calling with Mutect2 - Kristina's project"
__author__ = "Lars Andersen <larsmew@gmail.com> and Kristina Koldby"
__date__ = "16/12/2019"
__version__ = "1.0"

import time, os, sys, glob

#########################################################
####                       Input                     ####
#########################################################
# Matched tumour-normal samples information
configfile: "../samples-pairwise.yaml"
# configfile: "../samples-matched.yaml"

interval_vcf_dir = "/work/G65-2017-Kidstage/Connor/mutect2_somatic_variants/vcf/vcf_filterFlag_PASS/"
# force_intervals_suffix = "_mutect2_matched_filterFlag_norm_PASSonly_merged.vcf.gz"
force_intervals_suffix = "_somatic_mutect2_filterFlag_PASSonly.vcf.gz"

# Explicit paths for external input files
resource_path = "/work/sduvarcall/resources/hg38/"
ref = resource_path+"Homo_sapiens_assembly38.fasta"
gnomad = resource_path+"af-only-gnomad.hg38.vcf.gz"
interval_list = resource_path+"MedExome_hg38_capture_targets.interval_list"
target_regions = resource_path+"MedExome_target_regions/target_regions/"
common_variants = resource_path+"small_exac_common_3.hg38.vcf.gz"
# Panel of normals location
# pon_location = "somatic_panel_of_normals/"
pon_1000g = resource_path+"1000g_pon.hg38.vcf.gz"


SAMPLES = [sample for sample in config]
TUMORS = [config[sample]["tumor"] for sample in config]
NORMALS = [config[sample]["normal"] for sample in config]
# SAMPLES = ["ECV2-29"]
# TUMORS = ["ECV2-29-biopsi-merged_tumor_tagseq-medexome-deep-seq","ECV2-29-opA-merged_tumor_tagseq-medexome-deep-seq","ECV2-29-opB-merged_tumor_tagseq-medexome-deep-seq","ECV2-29-plasma171124_tumor_tagseq-medexome","ECV2-29-plasma180119_tumor_tagseq-medexome","ECV2-29-plasma180619_tumor_tagseq-medexome"]
# NORMALS = ["ECV2-29-blod_normal_tagseq-medexome"]

print(SAMPLES)
print(TUMORS)
print(NORMALS)

# Sample information
# PAIR, = glob_wildcards("{sample}_tumor.connor.recalibrated.bam")
# # PAIR = ["ECV2-4-plasma180205"]
# PAIR = [pair+"_tumor" for pair in PAIR]
# print(PAIR)

# INTERVALS, = glob_wildcards(target_regions+"{interval}.interval_list")
# INTERVALS = sorted(INTERVALS)
# print(INTERVALS)
# sys.exit()

# SAMPLES = "ECV2-29-blod_normal.connor.recalibrated_for_pon.vcf"

#########################################################
####                   Input/Output                  ####
#########################################################
log_file = "log_mutect2_joint_force_calling.txt"
output_somatic = "mutect2_force_calling_pairwise/"


#########################################################
####                       Setup                     ####
#########################################################
# Timing
totim = time.time()
timeFormat = "%Y_%m_%d:%X" # year, month, day, time H:M:S

# Memory
# mem = "-Xmx12g" # login nodes - be careful not running too many jobs at once!
mem = "-Xmx24g" # slim nodes
# mem = "-Xmx32g" # Fat nodes
# mem = "-Xmx64g"

onstart:
    # shell("mkdir -p "+output_somatic+"split/{{vcf,bam,f1r2}}_split")
    shell("mkdir -p "+output_somatic+"vcf_{{files,filterFlag,filterFlag_PASS,filterFlag_norm_PASS}}")
    # shell("mkdir -p "+output_somatic+"stats/")
    shell("mkdir -p "+output_somatic+"realigned_bam/")


#########################################################
####                  Run All Rules                  ####
#########################################################
'''
Rule all
'''
# print(expand(output_somatic+"{sample}_somatic_mutect2_filterFlag.vcf.gz", sample=SAMPLES))
# print(expand(output_somatic+"vcf_filterFlag_norm_decomposed/{sample}_somatic_mutect2_filterFlag_norm.vcf.gz", sample=SAMPLES))
# sys.exit()

rule all_pairs:
    input:
        # expand(output_somatic+"vcf_filterFlag_PASS/{sample}_somatic_mutect2_filterFlag_PASSonly.vcf.gz", sample=SAMPLES)
        expand(output_somatic+"vcf_final/{sample}_somatic_mutect2_filterFlag_norm_final.vcf.gz", sample=SAMPLES)
        # expand(output_somatic+"vcf_filterFlag_norm/{sample}_somatic_mutect2_filterFlag_norm.vcf.gz", sample=SAMPLES)
        # expand(output_somatic+"vcf_filterFlag_norm_PASS/{sample}_somatic_mutect2_filterFlag_norm_PASSonly.vcf.gz", sample=SAMPLES)
        # expand(output_somatic+"vcf_filterFlag_norm_decomposed/{sample}_somatic_mutect2_filterFlag_norm.vcf.gz", sample=SAMPLES)
        # expand(output_somatic+"vcf_files/{sample}_somatic_mutect2.vcf.gz", sample=SAMPLES),
        # expand(output_somatic+"{sample}_somatic_mutect2_filterFlag.vcf.gz", sample=SAMPLES),
        # expand(output_somatic+"{sample}_somatic_mutect2.bam", sample=SAMPLES)
        # expand(output_somatic+"{sample}_somatic_mutect2.vcf", sample=SAMPLES),
        # expand(output_somatic+"{sample}_merged_contamination.table", sample=SAMPLES)


##########################################################
####  Call Somatic Variants using Mutect2 on matched  ####
####   Tumor-Normal samples on per chromesome basis   ####
##########################################################
# Assumes PON already calculated in other mutect2 script.
'''
Mutect2 on matched Tumor-Normal samples
'''
rule Mutect2_joint_calling_forced:
    input:
        normal=lambda wildcards: expand("bam/{normal}.connor.recalibrated.bam", normal=config[wildcards.sample]["normal"]),
        tumors=lambda wildcards: expand("bam/{tumor}.connor.recalibrated.bam", tumor=config[wildcards.sample]["tumor"]),
        # pon=pon_location+"pon.vcf.gz",
        # intervals=target_regions+"{interval}.interval_list"
        # force_intervals=output_somatic+"../mutect2_somatic_variants/vcf/vcf_filterFlag_PASS/{sample}"+force_intervals_suffix
        # force_intervals=lambda wildcards: expand(output_somatic+"../mutect2_somatic_variants/vcf/vcf_filterFlag_PASS/{{sample}}_vs_{normal}"+force_intervals_suffix, normal=config[wildcards.sample]["normal"])
        force_intervals=lambda wildcards: expand(interval_vcf_dir+"{{sample}}_vs_{normal}"+force_intervals_suffix, normal=config[wildcards.sample]["normal"])
    output:
        vcf=output_somatic+"vcf_files/{sample}_somatic_mutect2.vcf.gz",
        idx=output_somatic+"vcf_files/{sample}_somatic_mutect2.vcf.gz.tbi",
        realigned_bam=output_somatic+"realigned_bam/{sample}_somatic_mutect2.bam",
        vcf_stats=output_somatic+"vcf_files/{sample}_somatic_mutect2.vcf.gz.stats",
        f1r2=output_somatic+"f1r2/{sample}_f1r2.tar.gz"
    params:
        tumors=lambda wildcards: expand("-I bam/{tumor}.connor.recalibrated.bam", tumor=config[wildcards.sample]["tumor"]),
        normal_name=lambda wildcards: expand("-normal {normal}", normal=config[wildcards.sample]["normal"])
    resources: cpus=4, mem_mb=24000
    shell:
        """
        gatk --java-options -Xmx{resources.mem_mb}g Mutect2 \
        -R {ref} \
        {params.tumors} \
        -I {input.normal} \
        {params.normal_name} \
        -pon {pon_1000g} \
        --germline-resource {gnomad} \
        --native-pair-hmm-threads {resources.cpus} \
        --f1r2-tar-gz {output.f1r2} \
        -L {input.force_intervals} \
        -alleles {input.force_intervals} \
        -bamout {output.realigned_bam} \
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
#         vcf_subfile=expand(output_somatic+"split/vcf_split/{{sample}}_somatic_mutect2__{interval}__split.vcf.gz", interval=INTERVALS)
#     output:
#         vcf=output_somatic+"vcf_files/{sample}_somatic_mutect2.vcf.gz",
#         idx=output_somatic+"vcf_files/{sample}_somatic_mutect2.vcf.gz.tbi",
#     params:
#         vcf_subfile=expand("-I "+output_somatic+"split/vcf_split/{{sample}}_somatic_mutect2__{interval}__split.vcf.gz", interval=INTERVALS)
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
#         bam_subfile=expand(output_somatic+"split/bam_split/{{sample}}__{interval}_mutect2.bam", interval=INTERVALS)
#   output:
#         bam=output_somatic+"realigned_bam/{sample}_somatic_mutect2.bam"
#   params:
#         bam_subfile=expand("-I "+output_somatic+"split/bam_split/{{sample}}__{interval}_mutect2.bam", interval=INTERVALS)
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
#         vcf_subfile=expand(output_somatic+"split/vcf_split/{{sample}}_somatic_mutect2__{interval}__split.vcf.gz.stats", interval=INTERVALS)
#     output:
#         vcf_stats=output_somatic+"vcf_files/{sample}_somatic_mutect2.vcf.gz.stats",
#     params:
#         vcf_subfile=expand("-stats "+output_somatic+"split/vcf_split/{{sample}}_somatic_mutect2__{interval}__split.vcf.gz.stats", interval=INTERVALS)
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
####       from formalin Formalin-Fixed Paraffin-    ####
####       Embedded (FFPE) samples - i.e. not needed ####
####       for frozen tissue                         ####
#########################################################
rule learnReadOrientationModel:
    input:
        f1r2=output_somatic+"f1r2/{sample}_f1r2.tar.gz"
    output:
        f1r2_model=output_somatic+"read_orientation_models/{sample}_read-orientation-model.tar.gz"
    shell:
        """
        gatk LearnReadOrientationModel \
        -I {input.f1r2} \
        -O {output.f1r2_model}
        """

#########################################################
####           Create Contamination table            ####
#########################################################
### Sæt sammen til en metode, som tager højde for med og uden normal... (Se gatk4 scripts)
rule GetPileupSummaries:
    input:
        bam="bam/{sample}.connor.recalibrated.bam"
    output:
        pileup=output_somatic+"pileups/{sample}_pileup.table"
    shell:
        """
        gatk --java-options {mem} GetPileupSummaries \
        -I {input.bam} \
        -V {common_variants} \
        -L {common_variants} \
        -O {output}
        """

rule CalculateContamination:
    input:
        normal=output_somatic+"pileups/{normal}_pileup.table",
        tumor=output_somatic+"pileups/{tumor}_pileup.table"
    output:
        contamination=output_somatic+"contaminations_tables/{tumor}_vs_{normal}_contamination.table"
    shell:
        """
        gatk --java-options {mem} CalculateContamination \
        -I {input.tumor} \
        -matched {input.normal} \
        -O {output}
        """

rule JoinContaminationTables:
    input:
        cont_tables=lambda wildcards: expand(output_somatic+"contaminations_tables/{tumor}_vs_{normal}_contamination.table", normal=config[wildcards.sample]["normal"], tumor=config[wildcards.sample]["tumor"])
    output:
        joined_table=output_somatic+"contaminations_tables/{sample}_merged_contamination.table"
    # params:
    #     sample_name=expand('{sample}', sample=SAMPLES)
    shell:
        """
        echo -e 'sample\tcontamination\terror' > {output};
        awk -F'\t' 'FNR == 2' {input} >> {output}
        """
        # awk -F'\t' 'FNR == 2' {params.sample_name}*contamination.table >> {output}


#########################################################
####              Filter Mutect2 Calls               ####
#########################################################
# rule FilterMutectCalls:
#     input:
#         vcf=output_somatic+"{sample}_somatic_mutect2.vcf",
#         idx=output_somatic+"{sample}_somatic_mutect2.vcf.idx",
#         stats=output_somatic+"{sample}_somatic_mutect2.vcf.stats",
#         vcf_stats=output_somatic+"{sample}_somatic_mutect2.vcf.stats",
#         contamination=output_somatic+"{sample}_merged_contamination.table",
#         read_orientation=output_somatic+"{sample}_read-orientation-model.tar.gz"
#     output:
#         vcf=output_somatic+"{sample}_somatic_mutect2_filtered.vcf",
#         idx=output_somatic+"{sample}_somatic_mutect2_filtered.vcf.idx"
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

rule FilterMutectCalls:
    input:
        vcf=output_somatic+"vcf_files/{sample}_somatic_mutect2.vcf.gz",
        idx=output_somatic+"vcf_files/{sample}_somatic_mutect2.vcf.gz.tbi",
        bam=output_somatic+"realigned_bam/{sample}_somatic_mutect2.bam",
        stats=output_somatic+"vcf_files/{sample}_somatic_mutect2.vcf.gz.stats",
        # contamination=output_somatic+"contaminations_tables/{sample}_merged_contamination.table",
        read_orientation=output_somatic+"read_orientation_models/{sample}_read-orientation-model.tar.gz"
    output:
        vcf=output_somatic+"vcf_filterFlag/{sample}_somatic_mutect2_filterFlag.vcf.gz",
        tbi=output_somatic+"vcf_filterFlag/{sample}_somatic_mutect2_filterFlag.vcf.gz.tbi",
        # vcf_pass=output_somatic+"vcf_filterFlag_PASS/{sample}_somatic_mutect2_filterFlag_PASSonly.vcf.gz",
        # tbi_pass=output_somatic+"vcf_filterFlag_PASS/{sample}_somatic_mutect2_filterFlag_PASSonly.vcf.gz.tbi"
    shell:
        """
        gatk --java-options {mem} FilterMutectCalls \
        -R {ref} \
        -V {input.vcf} \
        --orientation-bias-artifact-priors {input.read_orientation} \
        --stats {input.stats} \
        -L {interval_list} \
        -O {output.vcf}
        """
        # --contamination-table {input.contamination} \

rule vcf_pass_only:
    input:
        vcf=output_somatic+"vcf_filterFlag/{sample}_somatic_mutect2_filterFlag.vcf.gz",
    output:
        vcf_pass=output_somatic+"vcf_filterFlag_PASS/{sample}_somatic_mutect2_filterFlag_PASSonly.vcf.gz",
    shell:
        """
        bcftools view -i "%FILTER='PASS' | %FILTER='.'" {input.vcf} | bgzip -c > {output.vcf_pass}
        tabix -p vcf {output.vcf_pass}
        """

rule vcf_norm_decompose_bcftools:
    input:
        vcf=output_somatic+"vcf_filterFlag/{sample}_somatic_mutect2_filterFlag.vcf.gz"
    output:
        norm_vcf=output_somatic+"vcf_filterFlag_norm/{sample}_somatic_mutect2_filterFlag_norm.vcf.gz"
    shell:
        """
        bcftools norm -m-any {input.vcf} | \
            bcftools norm -Ou --check-ref w -f {ref} | \
            bcftools annotate -Oz -x 'ID' -I +'%CHROM:%POS:%REF:%ALT' > {output.norm_vcf}
        bcftools index -t {output.norm_vcf}
        """

rule vcf_norm_decompose_PASS_bcftools:
    input:
        vcf=output_somatic+"vcf_filterFlag/{sample}_somatic_mutect2_filterFlag.vcf.gz"
    output:
        norm_vcf=output_somatic+"vcf_filterFlag_norm_PASS/{sample}_somatic_mutect2_filterFlag_norm_PASSonly.vcf.gz"
    shell:
        """
        bcftools norm -m-any {input.vcf} | \
            bcftools norm -Ou --check-ref w -f {ref} | \
            bcftools filter -Ou -i 'FILTER="PASS"' | \
            bcftools annotate -Oz -x 'ID' -I +'%CHROM:%POS:%REF:%ALT' > {output.norm_vcf}
        bcftools index -t {output.norm_vcf}
        """

rule intersect_forcedVCF_mergedVCF:
    input:
        norm_vcf=output_somatic+"vcf_filterFlag_norm/{sample}_somatic_mutect2_filterFlag_norm.vcf.gz",
        # force_intervals=output_somatic+"union_vcf/{sample}"+force_intervals_suffix
        force_intervals=lambda wildcards: expand(interval_vcf_dir+"{{sample}}_vs_{normal}"+force_intervals_suffix, normal=config[wildcards.sample]["normal"])
    output:
        final_vcf=output_somatic+"vcf_final/{sample}_somatic_mutect2_filterFlag_norm_final.vcf.gz"
    shell:
        """
        bcftools isec -n+2 -w1 -Oz {input.norm_vcf} {input.force_intervals} > {output.final_vcf}
        bcftools index -t {output.final_vcf}
        """

# rule vcf_norm_decompose:
#     input:
#         vcf=output_somatic+"vcf_filterFlag/{sample}_somatic_mutect2_filterFlag.vcf.gz"
#     output:
#         norm_vcf=output_somatic+"vcf_filterFlag_norm_decomposed/{sample}_somatic_mutect2_filterFlag_norm.vcf.gz"
#     shell:
#         """
#         zless {input.vcf} \
#         | vt decompose -s - \
#         | vt normalize -r {ref} - \
#         | bgzip -c > {output.norm_vcf}
#         """
        # vcf index: tabix -p vcf {output.norm_vcf}
        # Alternative to vt
        # bcftools norm -m-any {input.vcf} > {output.vcf}

# rule vcf_norm_pass_only:
#     input:
#         norm_vcf=output_somatic+"vcf_filterFlag_norm_decomposed/{sample}_somatic_mutect2_filterFlag_norm.vcf.gz"
#     output:
#         vcf_pass=output_somatic+"vcf_filterFlag_norm_decomposed/{sample}_somatic_mutect2_filterFlag_norm_PASSonly.vcf.gz"
#     shell:
#         """
#         zcat {input.norm_vcf} | egrep '^#|PASS' | bgzip -c > {output.vcf_pass}
#         tabix -p vcf {output.vcf_pass}
#         """