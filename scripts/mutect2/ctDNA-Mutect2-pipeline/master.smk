__title__ = "Pipeline for Tumor Signal in Plasma using Mutect2"
__author__ = "Lars Andersen <larsmew@gmail.com>"
__date__ = "13/11/2021"
__version__ = "1.0"

import time, os, sys, glob

#########################################################
####                       Input                     ####
#########################################################
# Matched tumour-normal samples information
configfile: "/work/G65-2017-Kidstage/samples-perType.yaml"

# Explicit paths for external input files
resource_path = "/work/sduvarcall/resources/hg38/"
ref = resource_path+"Homo_sapiens_assembly38.fasta"
gnomad = resource_path+"af-only-gnomad.hg38.vcf.gz"
interval_list = resource_path+"MedExome_hg38_capture_targets.interval_list"
target_regions = resource_path+"MedExome_target_regions/target_regions/"
common_variants = resource_path+"small_exac_common_3.hg38.vcf.gz"
# Panel of normals location
pon_1000g = resource_path+"1000g_pon.hg38.vcf.gz"

#########################################################
####                      Output                     ####
#########################################################
main_output = "mutect2_ctDNA_pipeline/"
mutect2_matched_output = main_output+"mutect2_matched/"
mutect2_joint_calling_output = main_output+"mutect2_joint_calling/"
mutect2_force_calling_output = main_output+"mutect2_matched_force-calling/"


#########################################################
####                       Setup                     ####
#########################################################
# Timing
totim = time.time()
timeFormat = "%Y_%m_%d:%X" # year, month, day, time H:M:S

# Memory
mem = "-Xmx24g" # slim nodes
# mem = "-Xmx32g" # Fat nodes

# Samples
SAMPLES = [sample for sample in config]
TUMORS = [config[sample]["tumor"] for sample in config]
NORMALS = [config[sample]["normal"] for sample in config]
PLASMAS = [config[sample]["plasma"] for sample in config]
print(SAMPLES)

# Intervals - to split computations
INTERVALS, = glob_wildcards(target_regions+"{interval}.interval_list")
INTERVALS = sorted(INTERVALS)

# force_intervals_suffix = "_mutect2_matched_filterFlag_norm_PASSonly_merged.vcf.gz"
force_intervals_suffix = "_mutect2_matchedAndJoint_filterFlag_norm_PASSonly_merged.vcf.gz"

def get_plasma_pool(wildcards):
    plasmaPool = [sample for sample in SAMPLES if sample != wildcards.sample]
    print(plasmaPool)
    return plasmaPool

#########################################################
####                     Rule all                    ####
#########################################################
rule all:
    input:
        # expand("{output}vcf_filterFlag_PASS/{sample}_somatic_mutect2_filterFlag_PASSonly.vcf.gz", 
        #        sample=SAMPLES,
        #        output=mutect2_joint_calling_output)
        # expand("{output}mutect2_{sample}_finished.txt",
        #        sample=SAMPLES,
        #        output=mutect2_matched_output)
        # expand(main_output+"union_vcf/{sample}_mutect2_matchedAndJoint_filterFlag_norm_PASSonly_merged.vcf.gz", sample=SAMPLES)
        # expand(mutect2_force_calling_output+"vcf_final/{sample}_somatic_mutect2_filterFlag_norm_final.vcf.gz", sample=SAMPLES)
        [expand("{output}completed/mutect2_{sample}_plasmaPool_{sampleShort}_finished.txt",
               sample=sample,
               sampleShort=sample.split("_")[0],
            #    output=mutect2_force_calling_output) for sample in ["G74-C_thruplex-tag-seq-hv-nimblegen-medexome_HVFFJDSXY"]],
               output=mutect2_force_calling_output) for sample in SAMPLES],
        expand("{output}completed/mutect2_{sample}_finished.txt",
               sample=SAMPLES,
               output=mutect2_force_calling_output),
        # expand("{output}completed/mutect2_{sample}_tumorOnly_finished.txt",
        #        sample=SAMPLES,
        #        output=mutect2_force_calling_output),
        # expand("{output}completed/mutect2_{sample}_testAFcalc_finished.txt",
        #        sample=SAMPLES,
        #        output=mutect2_force_calling_output),
        # [expand(mutect2_force_calling_output+"vcf_final/maf_files/{tumor}_vs_{normal}__somatic_mutect2-forced_filterFlag_norm_final.maf",
        #     normal=config[fam]["normal"],
        #     tumor=config[fam]["tumor"]) for fam in PAIR],


rule compute_tumor_and_plasma:
    input:
        # lambda wildcards: expand(mutect2_force_calling_output+"vcf_final/{tumor}_vs_{normal}_somatic_mutect2-forced_filterFlag_norm_final.vcf.gz",
        lambda wildcards: expand(mutect2_force_calling_output+"vcf_final/maf_files/{tumor}_vs_{normal}__somatic_mutect2-forced_filterFlag_norm_final.maf",
               tumor=config[wildcards.sample]["tumor"],
               normal=config[wildcards.sample]["normal"]),
        lambda wildcards: expand(mutect2_force_calling_output+"vcf_final/maf_files/{plasma}_vs_{normal}__somatic_mutect2-forced_filterFlag_norm_final.maf",
               plasma=config[wildcards.sample]["plasma"],
               normal=config[wildcards.sample]["normal"]),
        # maf_dir=mutect2_force_calling_output+"vcf_final/maf_files/"
    output:
        # maf_out=main_output+
        "{output}completed/mutect2_{sample}_finished.txt"
    shell:
        """
        echo "SUCCESS" > {output}
        """

rule compute_tumor_and_plasma_tumorOnly:
    input:
        # lambda wildcards: expand(mutect2_force_calling_output+"vcf_final/{tumor}_vs_{normal}_somatic_mutect2-forced_filterFlag_norm_final.vcf.gz",
        lambda wildcards: expand(mutect2_force_calling_output+"vcf_final_tumorOnly/maf_files/{tumor}__somatic_mutect2-forced_filterFlag_norm_final_tumorOnly.maf",
               tumor=config[wildcards.sample]["tumor"]),
        lambda wildcards: expand(mutect2_force_calling_output+"vcf_final_tumorOnly/maf_files/{plasma}__somatic_mutect2-forced_filterFlag_norm_final_tumorOnly.maf",
               plasma=config[wildcards.sample]["plasma"]),
        # maf_dir=mutect2_force_calling_output+"vcf_final/maf_files/"
    output:
        # maf_out=main_output+
        "{output}completed/mutect2_{sample}_tumorOnly_finished.txt"
    shell:
        """
        echo "SUCCESS" > {output}
        """

# rule compute_tumor_and_plasma_testAFcalc:
#     input:
#         # lambda wildcards: expand(mutect2_force_calling_output+"vcf_final/{tumor}_vs_{normal}_somatic_mutect2-forced_filterFlag_norm_final.vcf.gz",
#         lambda wildcards: expand(mutect2_force_calling_output+"vcf_final_testAFcalc/maf_files/{tumor}_vs_{normal}__somatic_mutect2-forced_filterFlag_norm_final.maf",
#                tumor=config[wildcards.sample]["tumor"],
#                normal=config[wildcards.sample]["normal"]),
#         lambda wildcards: expand(mutect2_force_calling_output+"vcf_final_testAFcalc/maf_files/{plasma}_vs_{normal}__somatic_mutect2-forced_filterFlag_norm_final.maf",
#                plasma=config[wildcards.sample]["plasma"],
#                normal=config[wildcards.sample]["normal"]),
#         # maf_dir=mutect2_force_calling_output+"vcf_final/maf_files/"
#     output:
#         # maf_out=main_output+
#         "{output}completed/mutect2_{sample}_testAFcalc_finished.txt"
#     shell:
#         """
#         echo "SUCCESS" > {output}
#         """

rule compute_plasmaPool:
    input:
        # lambda wildcards: expand(mutect2_force_calling_output+"vcf_final/{tumor}_vs_{normal}_somatic_mutect2-forced_filterFlag_norm_final.vcf.gz",
        mafs=[expand(mutect2_force_calling_output+"plasmaPool/{{sample}}/vcf_final/plasmaPool_{{sampleShort}}_maf/{plasma}_vs_{normal}__somatic_mutect2-forced_filterFlag_norm_final.maf",
               plasma=config[s]["plasma"],
               normal=config[s]["normal"]) for s in SAMPLES],
        # maf_dir=directory(mutect2_force_calling_output+"plasmaPool/{sample}/vcf_final/maf_files/")
    output:
        # maf_out=directory(main_output+"maf_files/{sample}/plasmaPool/"),
        completed="{output}completed/mutect2_{sample}_plasmaPool_{sampleShort}_finished.txt"
    shell:
        """
        echo "SUCCESS" > {output.completed}
        """
        # cp -r {input.maf_dir} {output.maf_out}


#########################################################
####                       Rules                     ####
#########################################################

include: "rules/mutect2-joint-calling-gatk4-wes-somatic-hg38.smk"
include: "rules/mutect2-matched-gatk4-wes-somatic-hg38.smk"
include: "rules/mutect2-common-rules.smk"

include: "rules/Union-vcf-files.smk"

include: "rules/mutect2-matched-FORCE-calling-gatk4-wes-somatic-hg38.smk"
# include: "rules/mutect2-matched-FORCE-calling-gatk4-wes-somatic-hg38-plasmaPool.smk"
# include: "rules/mutect2-matched-FORCE-calling-gatk4-wes-somatic-hg38_tumorOnly.smk"
# include: "rules/mutect2-matched-FORCE-calling-gatk4-wes-somatic-hg38-testAFcalc.smk"
include: "rules/vcf2maf.smk"