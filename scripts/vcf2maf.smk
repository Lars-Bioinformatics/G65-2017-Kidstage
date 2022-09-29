__title__ = "vcf2maf"
__author__ = "Lars Andersen <larsmew@gmail.com>"
__date__ = "02/06/2021"
__version__ = "1.0"

import time

#########################################################
####                       Input                     ####
#########################################################
project = "/work/G65-2017-Kidstage/"
# Matched tumour-normal samples information
configfile: project+"samples-matched.yaml"
# configfile: project+"samples.yaml"
# configfile: "somatic_matched_samples_oneSample.yaml"

# Explicit paths for external input files
ref = "/work/sduvarcall/G37-2016-BRCAX-Genomes/reference_files_GRCh38/core_ref_GRCh38_hla_decoy_ebv/genome.fa"
# ref = "/work/sduvarcall/resources/hg38/Homo_sapiens_assembly38.fasta"
vep_data = "/work/sduvarcall/ensembl_vep/vep"
vep_path = "/work/miniconda3/envs/vep/bin"


# ref_build = "GRCh37"
ref_build = "GRCh38"


#########################################################
####                      Output                     ####
#########################################################
# output_mutect2 = "/work/Data/Connor/mutect2_vcf/"
output_mutect2 = project+"MarkDuplicates/mutect2_tumor_only_mode/vcf_filterFlag/"
# output_mutect2 = project+"WES/mutect2_somatic_variants/mutect2_vcf_filterFlag/"
# output_mutect2_pass = project+"Connor/Mutect2_somatic_variants/mutect2_vcf_filterFlag_pass/"

output_varscan2 = project+"WES/varscan_somatic_joint_calling/varscan_vcf_somatic_filtered/"
output_brass = project+"brass_ascat/"
# output_varscan2 = "/work/Data/Connor/varscan_somatic_joint_calling_q30/"

#########################################################
####               Sample Information                ####
#########################################################
# Mutect2 sample information
PAIR = [pair for pair in config]

# Varscan2 sample information
# VARSCAN_SAMPLES = glob_wildcards(output_varscan2+"{sample}_cns_varscan2.vcf.gz")
# print(VARSCAN_SAMPLES)

# All samples from Varscan2 config
# TUMORS = [x for y in [config[sample]["all"] for sample in config] for x in y]
# print(TUMORS)

#########################################################
####                       Setup                     ####
#########################################################
# Timing
totim = time.time()
timeFormat = "%Y_%m_%d:%X" # year, month, day, time H:M:S

# Memory
mem = "-Xmx12g"

#########################################################
####  Define workflow start, stop and error actions  ####
#########################################################
# onstart:
    # shell("mkdir -p "+output_mutect2)
    # shell("mkdir -p "+output_mutect2+"maf_files/")
    # shell("mkdir -p "+output_mutect2+"vcf_annotated/")
    # shell("mkdir -p "+output_mutect2+"mutect2_maf/")
    # shell("mkdir -p "+output_mutect2+"mutect2_maf_PASSonly/")
    # shell("mkdir -p "+output_varscan2+"vcf_annotated/")
    # shell("mkdir -p "+output_brass+"vcf_annotated/")
    # shell("mkdir -p "+output_varscan2+"varscan2_maf/")
    # shell("mkdir -p "+output_varscan2+"varscan2_maf_somatic_pon-filtered/")
    # shell("mkdir -p "+output_varscan2+"varscan2_maf_somatic_1000g-pon-filtered/")
    # shell("mkdir -p "+output_varscan2+"varscan2_maf_somatic_1000g-pon_gnomad-filtered/")


#########################################################
####                  Run All Rules                  ####
#########################################################
'''
Rule all
'''
rule all_pairs:
    input:
        # Tumor-only vcf
        # expand(output_mutect2+"maf_files/{tumor}__{extra}_norm_tumorOnly.maf", tumor=TUMORS, extra="somatic_mutect2_filterFlag"),
        # Multi-sample vcf
        # [expand("mutect2_maf/{tumor}_vs_{normal}__somatic_mutect2_filterFlag_fromMultiSampleVcf.maf",
        #     normal=config[fam]["normal"],
        #     tumor=config[fam]["tumor"]) for fam in PAIR],
        # [expand("mutect2_maf/{tumor}_vs_{normal}__somatic_mutect2_filterFlag_norm_final_fromMultiSampleVcf.maf",
        #     normal=config[fam]["normal"],
        #     tumor=config[fam]["tumor"]) for fam in PAIR],
        # [expand(output_mutect2+"maf_files/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm.maf",
        #     normal=config[fam]["normal"],
        #     tumor=config[fam]["tumor"]) for fam in PAIR],
        # [expand(output_brass+"maf_files/{tumor}_vs_{normal}.annot_norm.maf",
        # Multi-sample vcf - with plasma pool
        [expand("maf_files/{tumor}_plasmaPool/{sample}__somatic_mutect2_filterFlag_norm_final_fromMultiSampleVcf.maf",
            tumor=tumor,
            sample=config[tumor]["tumor"]) for tumor in config],
        #     normal=config[fam]["normal"],
        #     tumor=config[fam]["tumor"]) for fam in PAIR],
        # [expand(output_mutect2+"mutect2_maf_PASSonly/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm_PASSonly.maf",
        #     normal=config[fam]["normal"],
        #     tumor=config[fam]["tumor"]) for fam in PAIR],
        # [expand(output_varscan2+"varscan2_maf_somatic_1000g-pon-filtered/{tumor}_vs_{normal}_cns_varscan2_somatic_1000g-pon-filtered.maf", 
        #     normal=config[fam]["normal"],
        #     tumor=config[fam]["tumor"]) for fam in PAIR],
        # [expand(output_varscan2+"varscan2_maf_somatic_pon-filtered/{tumor}_cns_varscan2_somatic_pon-filtered.maf", 
        #     tumor=config[fam]["tumor"]) for fam in PAIR],
        # [expand(output_varscan2+"varscan2_maf_somatic_1000g-pon_gnomad-filtered/{tumor}_vs_{normal}_cns_varscan2_somatic_1000g-pon_gnomad-filtered.maf", 
        #     normal=config[fam]["normal"],
        #     tumor=config[fam]["tumor"]) for fam in PAIR],
        # [expand(output_varscan2+"varscan2_maf/{tumor}_vs_{normal}_cns_varscan2.maf",
        #     normal=config[fam]["normal"],
        #     tumor=config[fam]["tumor"]) for fam in PAIR],


#########################################################
####                  Run All Rules                  ####
#########################################################
rule vcf2maf_pairwise:
    input:
        vcf="{path}{tumor}_vs_{normal}_{extra}.vcf.gz"
    output:
        maf="{path}maf_files/{tumor}_vs_{normal}_{extra}_norm.maf",
        vep_vcf="{path}vcf_annotated/{tumor}_vs_{normal}_{extra}_norm.vep.vcf",
        vcf=temp("{path}vcf_annotated/{tumor}_vs_{normal}_{extra}_norm.vcf")
    params:
        new_tum_id = lambda wildcards: wildcards.tumor.split("_")[0],
        new_norm_id = lambda wildcards: str(wildcards.normal).split("_")[0]
    resources: cpus=4, mem_mb=24000 # Since internal command of vcf2maf use fork = 4
    shell:
        """
        bgzip -d -c {input.vcf} > {output.vcf}
        conda run --name vep \
        vcf2maf.pl \
        --input-vcf {output.vcf} \
        --output-maf {output.maf} \
        --vep-data {vep_data} \
        --vep-path {vep_path} \
        --ncbi-build {ref_build} \
        --ref-fasta {ref} \
        --tumor-id {params.new_tum_id} \
        --normal-id {params.new_norm_id} \
        --vcf-tumor-id {wildcards.tumor} \
        --vcf-normal-id {wildcards.normal}
        """

rule vcf2maf_multiSample:
    input:
        vcf=lambda wildcards: expand("{sample}_{{suffix}}.vcf.gz", sample=wildcards.normal[:5]+wildcards.normal[7:])
    output:
        maf="mutect2_maf/{tumor}_vs_{normal}__{suffix}_fromMultiSampleVcf.maf",
        vep_vcf="vcf_annotated/{tumor}_vs_{normal}__{suffix}_fromMultiSampleVcf.vep.vcf",
        vcf=temp("vcf_annotated/{tumor}_vs_{normal}__{suffix}_fromMultiSampleVcf.vcf")
    params:
        new_tum_id = lambda wildcards: wildcards.tumor.split("_")[0],
        new_norm_id = lambda wildcards: str(wildcards.normal).split("_")[0]
    resources: cpus=4, mem_mb=24000 # Since internal command of vcf2maf use fork = 4
    shell:
        """
        bgzip -d -c {input.vcf} > {output.vcf}
        conda run --name vep \
        vcf2maf.pl \
        --input-vcf {output.vcf} \
        --output-maf {output.maf} \
        --vep-data {vep_data} \
        --vep-path {vep_path} \
        --ncbi-build {ref_build} \
        --ref-fasta {ref} \
        --tumor-id {params.new_tum_id} \
        --normal-id {params.new_norm_id} \
        --vcf-tumor-id {wildcards.tumor} \
        --vcf-normal-id {wildcards.normal}
        """

rule vcf2maf_tumorOnly:
    input:
        vcf="{path}{tumor}_{extra}.vcf.gz"
    output:
        maf="{path}maf_files/{tumor}__{extra}_norm_tumorOnly.maf",
        vep_vcf="{path}vcf_annotated/{tumor}__{extra}_norm_tumorOnly.vep.vcf",
        vcf=temp("{path}vcf_annotated/{tumor}__{extra}_norm_tumorOnly.vcf")
    params:
        new_tum_id = lambda wildcards: wildcards.tumor.split("_")[0],
    resources: cpus=4, mem_mb=24000 # Since internal command of vcf2maf use fork = 4
    shell:
        """
        bgzip -d -c {input.vcf} > {output.vcf}
        conda run --name vep \
        vcf2maf.pl \
        --input-vcf {output.vcf} \
        --output-maf {output.maf} \
        --vep-data {vep_data} \
        --vep-path {vep_path} \
        --ncbi-build {ref_build} \
        --ref-fasta {ref} \
        --tumor-id {params.new_tum_id} \
        --vcf-tumor-id {wildcards.tumor}
        """

rule vcf2maf_tumor_plasmaPool:
    input:
        vcf="{tumor}_{suffix}.vcf.gz"
    output:
        maf="maf_files/{tumor}_plasmaPool/{sample}__{suffix}_fromMultiSampleVcf.maf",
        vep_vcf="vcf_annotated/{tumor}_plasmaPool/{sample}__{suffix}_fromMultiSampleVcf.vep.vcf",
        vcf=temp("vcf_annotated/{tumor}_plasmaPool/{sample}__{suffix}_fromMultiSampleVcf.vcf")
    params:
        new_tum_id = lambda wildcards: wildcards.sample.split("_")[0],
        normal = lambda wildcards: config[wildcards.tumor]["normal"],
        new_norm_id = lambda wildcards: config[wildcards.tumor]["normal"].split("_")[0]
    resources: cpus=2, mem_mb=12000 # Since internal command of vcf2maf use fork = 4
    shell:
        """
        bgzip -d -c {input.vcf} > {output.vcf}
        conda run --name vep \
        vcf2maf.pl \
        --input-vcf {output.vcf} \
        --output-maf {output.maf} \
        --vep-data {vep_data} \
        --vep-path {vep_path} \
        --ncbi-build {ref_build} \
        --ref-fasta {ref} \
        --tumor-id {params.new_tum_id} \
        --normal-id {params.new_norm_id} \
        --vcf-tumor-id {wildcards.sample} \
        --vcf-normal-id {params.normal}
        """

# rule vcf2maf_mutect2:
#     input:
#         vcf=output_mutect2+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag.vcf.gz"
#     output:
#         maf=output_mutect2+"mutect2_maf/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm.maf",
#         vep_vcf=output_mutect2+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm.vep.vcf",
#         vcf=temp(output_mutect2+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm.vcf")
#     #params:
#     #    vcf=temp(output_mutect2+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag.vcf")
#     resources: cpus=4, mem=24000 # Since internal command of vcf2maf use fork = 4
#     shell:
#         """
#         bgzip -d -c {input.vcf} > {output.vcf}
#         vcf2maf.pl \
#         --input-vcf {output.vcf} \
#         --output-maf {output.maf} \
#         --vep-data {vep_data} \
#         --vep-path {vep_path} \
#         --ncbi-build {ref_build} \
#         --ref-fasta {ref} \
#         --tumor-id {wildcards.tumor} \
#         --normal-id {wildcards.normal} \
#         --vcf-tumor-id {wildcards.tumor} \
#         --vcf-normal-id {wildcards.normal}
#         """

# rule vcf2maf_mutect2_multiSample:
#     input:
#         vcf=lambda wildcards: expand(output_mutect2+"{sample}_somatic_mutect2_filterFlag.vcf.gz", sample=wildcards.normal[:5]+wildcards.normal[7:])
#     output:
#         maf=output_mutect2+"mutect2_maf/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm_fromMultiSampleVcf.maf",
#         vep_vcf=output_mutect2+"vcf_annotated/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm_fromMultiSampleVcf.vep.vcf",
#         vcf=temp(output_mutect2+"vcf_annotated/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm_fromMultiSampleVcf.vcf")
#     params:
#         new_tum_id = lambda wildcards: wildcards.tumor.split("_")[0],
#         new_norm_id = lambda wildcards: str(wildcards.normal).split("_")[0]
#     resources: cpus=4, mem=24000 # Since internal command of vcf2maf use fork = 4
#     shell:
#         """
#         bgzip -d -c {input.vcf} > {output.vcf}
#         vcf2maf.pl \
#         --input-vcf {output.vcf} \
#         --output-maf {output.maf} \
#         --vep-data {vep_data} \
#         --vep-path {vep_path} \
#         --ncbi-build {ref_build} \
#         --ref-fasta {ref} \
#         --tumor-id {params.new_tum_id} \
#         --normal-id {params.new_norm_id} \
#         --vcf-tumor-id {wildcards.tumor} \
#         --vcf-normal-id {wildcards.normal}
#         """



# rule vcf2maf_varscan2:
#     input:
#         vcf=lambda wildcards: expand(output_varscan2+"{sample}_cns_varscan2.vcf.gz", sample=wildcards.normal[:5]+wildcards.normal[7:])
#     output:
#         maf=output_varscan2+"varscan2_maf/{tumor}_vs_{normal}_cns_varscan2.maf",
#         vep_vcf=output_varscan2+"vcf_annotated/{tumor}_vs_{normal}_cns_varscan2.vep.vcf",
#         vcf=temp(output_varscan2+"vcf_annotated/{tumor}_vs_{normal}_cns_varscan2.vcf")
#     params:
#         new_tum_id = lambda wildcards: wildcards.tumor.split("_")[0],
#         new_norm_id = lambda wildcards: str(wildcards.normal).split("_")[0]
#     resources: cpus=4, mem=24000 # Since internal command of vcf2maf use fork = 4
#     shell:
#         """
#         bgzip -d -c {input.vcf} > {output.vcf}
#         vcf2maf.pl \
#         --input-vcf {output.vcf} \
#         --output-maf {output.maf} \
#         --vep-data {vep_data} \
#         --vep-path {vep_path} \
#         --ncbi-build {ref_build} \
#         --ref-fasta {ref} \
#         --tumor-id {params.new_tum_id} \
#         --normal-id {params.new_norm_id} \
#         --vcf-tumor-id {wildcards.tumor} \
#         --vcf-normal-id {wildcards.normal}
#         """

# rule vcf2maf_varscan2_somtic_pon:
#     input:
#         vcf=output_varscan2+"{tumor}_cns_varscan2_somatic_pon-filtered.vcf.gz"
#     output:
#         maf=output_varscan2+"varscan2_maf_somatic_pon-filtered/{tumor}_cns_varscan2_somatic_pon-filtered.maf",
#         vep_vcf=output_varscan2+"vcf_annotated/{tumor}_cns_varscan2_somatic_pon-filtered.vep.vcf",
#         vcf=temp(output_varscan2+"vcf_annotated/{tumor}_cns_varscan2_somatic_pon-filtered.vcf")
#     params:
#         normal = lambda wildcards: config[wildcards.tumor]["normal"],
#         new_tum_id = lambda wildcards: wildcards.tumor.split("_")[0],
#         new_norm_id = lambda wildcards: config[wildcards.tumor]["normal"].split("_")[0]
#     resources: cpus=4, mem=24000 # Since internal command of vcf2maf use fork = 4
#     shell:
#         """
#         bgzip -d -c {input.vcf} > {output.vcf}
#         vcf2maf.pl \
#         --input-vcf {output.vcf} \
#         --output-maf {output.maf} \
#         --vep-data {vep_data} \
#         --vep-path {vep_path} \
#         --ncbi-build {ref_build} \
#         --ref-fasta {ref} \
#         --tumor-id {params.new_tum_id} \
#         --normal-id {params.new_norm_id} \
#         --vcf-tumor-id {wildcards.tumor} \
#         --vcf-normal-id {params.normal}
#         """


# rule vcf2maf_varscan2_somtic_1000g_pon:
#     input:
#         vcf=lambda wildcards: expand(output_varscan2+"{sample}_cns_varscan2_somatic_1000g-pon-filtered.vcf.gz", sample=wildcards.normal[:5]+wildcards.normal[7:])
#     output:
#         maf=output_varscan2+"varscan2_maf_somatic_1000g-pon-filtered/{tumor}_vs_{normal}_cns_varscan2_somatic_1000g-pon-filtered.maf",
#         vep_vcf=output_varscan2+"vcf_annotated/{tumor}_vs_{normal}_cns_varscan2_somatic_1000g-pon-filtered.vep.vcf",
#         vcf=temp(output_varscan2+"vcf_annotated/{tumor}_vs_{normal}_cns_varscan2_somatic_1000g-pon-filtered.vcf")
#     params:
#         new_tum_id = lambda wildcards: wildcards.tumor.split("_")[0],
#         new_norm_id = lambda wildcards: str(wildcards.normal).split("_")[0]
#     resources: cpus=4, mem=24000 # Since internal command of vcf2maf use fork = 4
#     shell:
#         """
#         bgzip -d -c {input.vcf} > {output.vcf}
#         vcf2maf.pl \
#         --input-vcf {output.vcf} \
#         --output-maf {output.maf} \
#         --vep-data {vep_data} \
#         --vep-path {vep_path} \
#         --ncbi-build {ref_build} \
#         --ref-fasta {ref} \
#         --tumor-id {params.new_tum_id} \
#         --normal-id {params.new_norm_id} \
#         --vcf-tumor-id {wildcards.tumor} \
#         --vcf-normal-id {wildcards.normal}
#         """

# rule vcf2maf_varscan2_somtic_pon_gnomad:
#     input:
#         vcf=lambda wildcards: expand(output_varscan2+"{sample}_cns_varscan2_somatic_1000g-pon_gnomad-filtered.vcf.gz", sample=wildcards.normal[:5]+wildcards.normal[7:])
#     output:
#         maf=output_varscan2+"varscan2_maf_somatic_1000g-pon_gnomad-filtered/{tumor}_vs_{normal}_cns_varscan2_somatic_1000g-pon_gnomad-filtered.maf",
#         vep_vcf=output_varscan2+"vcf_annotated/{tumor}_vs_{normal}_cns_varscan2_somatic_1000g-pon_gnomad-filtered.vep.vcf",
#         vcf=temp(output_varscan2+"vcf_annotated/{tumor}_vs_{normal}_cns_varscan2_somatic_1000g-pon_gnomad-filtered.vcf")
#     params:
#         new_tum_id = lambda wildcards: wildcards.tumor.split("_")[0],
#         new_norm_id = lambda wildcards: str(wildcards.normal).split("_")[0]
#     resources: cpus=4, mem=24000 # Since internal command of vcf2maf use fork = 4
#     shell:
#         """
#         bgzip -d -c {input.vcf} > {output.vcf}
#         vcf2maf.pl \
#         --input-vcf {output.vcf} \
#         --output-maf {output.maf} \
#         --vep-data {vep_data} \
#         --vep-path {vep_path} \
#         --ncbi-build {ref_build} \
#         --ref-fasta {ref} \
#         --tumor-id {params.new_tum_id} \
#         --normal-id {params.new_norm_id} \
#         --vcf-tumor-id {wildcards.tumor} \
#         --vcf-normal-id {wildcards.normal}
#         """


