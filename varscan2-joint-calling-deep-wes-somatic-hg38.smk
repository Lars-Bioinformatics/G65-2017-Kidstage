__title__ = "Pipeline for Somatic JOINT Variant Calling with Varscan2 - Kristina's project"
__author__ = "Lars Andersen <larsmew@gmail.com> and Kristina Koldby"
__date__ = "8/1/2020"
__version__ = "1.0"

import time, os, sys, glob

#########################################################
####                       Input                     ####
#########################################################
# Matched tumour-normal samples information
configfile: "../samples.yaml"
# configfile: "somatic_matched_samples_oneSample.yaml"

# Explicit paths for external input files
ref = "/work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.fasta"
target_regions = "/work/sduvarcall/Resources/hg38/MedExome_target_regions/target_regions/bed_files/"


# PATIENTS = ["ECV2-4", "ECV2-8", "ECV2-29", "ECV2-31", "ECV2-35", "PC1-10", "PC1-14", "PC1-18", "PON-all"]
# PATIENTS = ["G65-T28A-42A6_nimblegen-medexome_HYLKFDSXX"]
PATIENTS = [sample for sample in config if sample.startswith("G65")]
PATIENTS.append("PON-all")


INTERVALS, = glob_wildcards(target_regions+"{interval}.bed")
INTERVALS = sorted(INTERVALS)
# print(INTERVALS)
print (PATIENTS)


#########################################################
####                      Output                     ####
#########################################################
output_somatic = "varscan_somatic_joint_calling/"


#########################################################
####                       Setup                     ####
#########################################################
# Timing
totim = time.time()
timeFormat = "%Y_%m_%d:%X" # year, month, day, time H:M:S

# Memory
mem = "-Xmx12g" # login nodes
# mem = "-Xmx24g" # slim nodes
# mem = "-Xmx32g" # Fat nodes
# mem = "-Xmx64g"

onstart:
    shell("mkdir -p "+output_somatic)
    shell("mkdir -p "+output_somatic+"split")
    shell("mkdir -p "+output_somatic+"split_vcf")


#########################################################
####                  Run All Rules                  ####
#########################################################
'''
Rule all
'''
print(expand(output_somatic+"{sample}_cns_varscan2.vcf.gz", sample=PATIENTS))
# sys.exit()

rule all_pairs:
    input:
        expand(output_somatic+"{sample}_cns_varscan2.vcf.gz", sample=PATIENTS),
        expand(output_somatic+"{sample}_cns_varscan2_somatic_pon-filtered.vcf.gz", sample=PATIENTS),
        expand(output_somatic+"{sample}_cns_varscan2_somatic_pon-filtered_strict.vcf.gz", sample=PATIENTS),
        expand(output_somatic+"{sample}_cns_varscan2_somatic.vcf.gz", sample=PATIENTS),
        expand(output_somatic+"{sample}_somaticCall_varscan2_somatic_pon-filtered_snp.vcf.gz", sample="G65-T2-C99C_nimblegen-medexome_HYLKFDSXX"),
        expand(output_somatic+"{sample}_somaticCall_varscan2_somatic_pon-filtered_strict_snp.vcf.gz", sample="G65-T2-C99C_nimblegen-medexome_HYLKFDSXX"),
        expand(output_somatic+"{sample}_somaticCall_varscan2_somatic_pon-filtered_indel.vcf.gz", sample="G65-T2-C99C_nimblegen-medexome_HYLKFDSXX"),
        expand(output_somatic+"{sample}_somaticCall_varscan2_somatic_pon-filtered_strict_indel.vcf.gz", sample="G65-T2-C99C_nimblegen-medexome_HYLKFDSXX"),
        # expand(output_somatic+"{sample}_cns_varscan2_somatic_NormADzero.vcf.gz", sample=PATIENTS),
        # expand(output_somatic+"{sample}_sampleList.txt", sample=PATIENTS)


##########################################################
####  Create list of sample names                     ####
##########################################################
rule sampleList:
    output:
        output_somatic+"{sample}_sampleList.txt"
    params:
        sample_names = lambda wildcards: config[wildcards.sample]["all"]
    run:
        # print("_tagseq-medexome\n".join(SAMPLE_NAMES))
        # print(output[0])
        with open(output[0],"w") as f:
            f.write("\n".join([s for s in params[0]]))
            # f.write("\n".join([s+"_tagseq-medexome" for s in params[0]]))


##########################################################
####  Samtools mpileup                                ####
##########################################################
rule mpileup:
    input:
        bam = lambda wildcards: ["bam/"+b+".connor.recalibrated.bam" for b in config[wildcards.sample]["all"]],
        intervals=target_regions+"{interval}.bed"
    output:
        mpileups=output_somatic+"split/{sample}__{interval}__split.mpileup"
    shell:
        """
        samtools mpileup \
        -f {ref} \
        -l {input.intervals}\
        -d 1000000 \
        {input.bam} > {output}
        """
        # -B \


##########################################################
####  Call Somatic Variants using Varscan2 on matched ####
####   Tumor-Normal samples                           ####
##########################################################
'''
Multi-sample calling using Varscan2
'''
rule Varscan2:
    input:
        mpileups=output_somatic+"split/{sample}__{interval}__split.mpileup",
        sample_list=output_somatic+"{sample}_sampleList.txt"
    output:
        vcf=output_somatic+"split_vcf/{sample}_cns__{interval}__split_varscan2.vcf.gz",
        vcf_index=output_somatic+"split_vcf/{sample}_cns__{interval}__split_varscan2.vcf.gz.tbi"
    shell:
        """
        varscan mpileup2cns {input.mpileups} \
        --min-coverage 1 \
        --min-reads2 2 \
        --min-var-freq 0 \
        --output-vcf 1 \
        --min-avg-qual 20 \
        --variants 1 \
        --strand-filter 0 \
        --p-value 0.01 \
        --vcf-sample-list {input.sample_list} \
        | bgzip -c > {output.vcf}
        bcftools index {output.vcf} -o {output.vcf_index}
        """
        # Default p-value is 0.01


# rule index_vcf:
#     input:
#         vcf_subfile=output_somatic+"split_vcf/{sample}_cns__{interval}__split_varscan2.vcf"
#     output:
#         vcf_subfile_compressed=output_somatic+"split_vcf/{sample}_cns__{interval}__split_varscan2.vcf.gz",
#         vcf_index=output_somatic+"split_vcf/{sample}_cns__{interval}__split_varscan2.vcf.gz.tbi"
#     shell:
#         """
#         bgzip -c {input} > {output.vcf_subfile_compressed}
#         bcftools index {output.vcf_subfile_compressed} -o {output.vcf_index}
#         """

rule merge_Varscan2:
    input:
        vcf_subfile=expand(output_somatic+"split_vcf/{{sample}}_cns__{interval}__split_varscan2.vcf.gz", interval=INTERVALS),
        vcf_index=expand(output_somatic+"split_vcf/{{sample}}_cns__{interval}__split_varscan2.vcf.gz.tbi", interval=INTERVALS)
    output:
        vcf=output_somatic+"{sample}_cns_varscan2.vcf.gz",
        vcf_index=output_somatic+"{sample}_cns_varscan2.vcf.gz.tbi"
    params:
        temp=output_somatic+"{sample}_cns_varscan2.vcf"
    shell:
        """
        bcftools concat \
        {input.vcf_subfile} \
        --no-version \
        -o {params.temp};
        bgzip {params.temp};
        tabix -p vcf {output.vcf}
        """

        # """
        # bcftools concat \
        # {input.vcf_subfile} \
        # --no-version \
        # -o {output.vcf};
        # bgzip -c {output.vcf} > {output.vcf}.gz;
        # tabix -p vcf {output.vcf}.gz
        # """

### Remove variants found with more than one read i blood sample ###
### IMPORTANT: 'FORMAT/AD[0]' requires that blood sample is the first sample in joint vcf
### This is the case when blood sample is listed as first sample in category "all" i config file

rule somatic_filter:
    input:
        output_somatic+"{sample}_cns_varscan2.vcf.gz"
    output:
        output_somatic+"{sample}_cns_varscan2_somatic.vcf.gz"
    shell:
        """
        bcftools filter \
        --include 'FORMAT/AD[0] <= 1' \
        {input} \
        -o {output}
        """

rule somatic_filter_strict:
    input:
        output_somatic+"{sample}_cns_varscan2.vcf.gz"
    output:
        output_somatic+"{sample}_cns_varscan2_somatic_NormADzero.vcf.gz"
    shell:
        """
        bcftools filter \
        --include 'FORMAT/AD[0] = 0' \
        {input} \
        -o {output}
        """

rule pon_filtering:
    input:
        vcf = output_somatic+"{sample}_cns_varscan2.vcf.gz",
        pon = output_somatic+"PON-all_cns_varscan2.vcf.gz"
    output:
        output_somatic+"{sample}_cns_varscan2_somatic_pon-filtered.vcf.gz"
    shell:
        """
        bcftools isec \
        -C \
        -c none \
        -w1 \
        {input.vcf} \
        {input.pon} | \
        bcftools filter \
        --include 'FORMAT/AD[0] <= 1' \
        -o {output}
        """

rule pon_filtering_strict:
    input:
        vcf = output_somatic+"{sample}_cns_varscan2.vcf.gz",
        pon = output_somatic+"PON-all_cns_varscan2.vcf.gz"
    output:
        output_somatic+"{sample}_cns_varscan2_somatic_pon-filtered_strict.vcf.gz"
    shell:
        """
        bcftools isec \
        -C \
        -c none \
        -w1 \
        {input.vcf} \
        {input.pon} | \
        bcftools filter \
        --include 'FORMAT/AD[0] = 0' \
        -o {output}
        """






# VARSCAN SOMATIC MODE

##########################################################
####  Call Somatic Variants using Varscan2 somatic    ####
####  option on matched Tumor-Normal samples          ####
##########################################################
'''
Multi-sample calling using Varscan2
'''
rule Varscan2_somatic:
    input:
        mpileups=output_somatic+"split/{sample}__{interval}__split.mpileup",
        sample_list=output_somatic+"{sample}_sampleList.txt"
    output:
        vcf_snp=output_somatic+"split_vcf/{sample}_somatic__{interval}__split_varscan2_snp.vcf.gz",
        vcf_indel=output_somatic+"split_vcf/{sample}_somatic__{interval}__split_varscan2_indel.vcf.gz",
        # vcf_index=output_somatic+"split_vcf/{sample}_somatic__{interval}__split_varscan2.vcf.gz.tbi"
    shell:
        """
        varscan somatic {input.mpileups} \
        --mpileup 1 \
        --min-coverage 1 \
        --min-var-freq 0 \
        --output-vcf 1 \
        --strand-filter 0 \
        --p-value 0.01 \
        --output-snp {output.vcf_snp} \
        --output-indel {output.vcf_indel}
        """
        # Default p-value is 0.01


# rule index_vcf:
#     input:
#         vcf_subfile=output_somatic+"split_vcf/{sample}_cns__{interval}__split_varscan2.vcf"
#     output:
#         vcf_subfile_compressed=output_somatic+"split_vcf/{sample}_cns__{interval}__split_varscan2.vcf.gz",
#         vcf_index=output_somatic+"split_vcf/{sample}_cns__{interval}__split_varscan2.vcf.gz.tbi"
#     shell:
#         """
#         bgzip -c {input} > {output.vcf_subfile_compressed}
#         bcftools index {output.vcf_subfile_compressed} -o {output.vcf_index}
#         """

rule merge_Varscan2_somatic_snp:
    input:
        vcf_subfile=expand(output_somatic+"split_vcf/{{sample}}_somatic__{interval}__split_varscan2_snp.vcf.gz", interval=INTERVALS),
        # vcf_index=expand(output_somatic+"split_vcf/{{sample}}_somatic__{interval}__split_varscan2.vcf.gz.tbi", interval=INTERVALS)
    output:
        vcf=output_somatic+"{sample}_somatic_varscan2_snp.vcf.gz",
        vcf_index=output_somatic+"{sample}_somatic_varscan2_snp.vcf.gz.tbi"
    params:
        temp=output_somatic+"{sample}_somatic_varscan2_snp.vcf"
    shell:
        """
        bcftools concat \
        {input.vcf_subfile} \
        --no-version \
        -o {params.temp};
        bgzip {params.temp};
        tabix -p vcf {output.vcf}
        """

rule merge_Varscan2_somatic_indel:
    input:
        vcf_subfile=expand(output_somatic+"split_vcf/{{sample}}_somatic__{interval}__split_varscan2_indel.vcf.gz", interval=INTERVALS),
        # vcf_index=expand(output_somatic+"split_vcf/{{sample}}_somatic__{interval}__split_varscan2.vcf.gz.tbi", interval=INTERVALS)
    output:
        vcf=output_somatic+"{sample}_somatic_varscan2_indel.vcf.gz",
        vcf_index=output_somatic+"{sample}_somatic_varscan2_indel.vcf.gz.tbi"
    params:
        temp=output_somatic+"{sample}_somatic_varscan2_indel.vcf"
    shell:
        """
        bcftools concat \
        {input.vcf_subfile} \
        --no-version \
        -o {params.temp};
        bgzip {params.temp};
        tabix -p vcf {output.vcf}
        """


rule pon_filtering_somaticCall_snp:
    input:
        vcf = output_somatic+"{sample}_somatic_varscan2_snp.vcf.gz",
        pon = output_somatic+"PON-all_cns_varscan2.vcf.gz"
    output:
        output_somatic+"{sample}_somaticCall_varscan2_somatic_pon-filtered_snp.vcf.gz"
    shell:
        """
        bcftools isec \
        -C \
        -c none \
        -w1 \
        {input.vcf} \
        {input.pon} | \
        bcftools filter \
        --include 'FORMAT/AD[0] <= 1' \
        -o {output}
        """

rule pon_filtering_strict_somaticCall_snp:
    input:
        vcf = output_somatic+"{sample}_somatic_varscan2_snp.vcf.gz",
        pon = output_somatic+"PON-all_cns_varscan2.vcf.gz"
    output:
        output_somatic+"{sample}_somaticCall_varscan2_somatic_pon-filtered_strict_snp.vcf.gz"
    shell:
        """
        bcftools isec \
        -C \
        -c none \
        -w1 \
        {input.vcf} \
        {input.pon} | \
        bcftools filter \
        --include 'FORMAT/AD[0] = 0' \
        -o {output}
        """ 

rule pon_filtering_somaticCall_indel:
    input:
        vcf = output_somatic+"{sample}_somatic_varscan2_indel.vcf.gz",
        pon = output_somatic+"PON-all_cns_varscan2.vcf.gz"
    output:
        output_somatic+"{sample}_somaticCall_varscan2_somatic_pon-filtered_indel.vcf.gz"
    shell:
        """
        bcftools isec \
        -C \
        -c none \
        -w1 \
        {input.vcf} \
        {input.pon} | \
        bcftools filter \
        --include 'FORMAT/AD[0] <= 1' \
        -o {output}
        """

rule pon_filtering_strict_somaticCall_indel:
    input:
        vcf = output_somatic+"{sample}_somatic_varscan2_indel.vcf.gz",
        pon = output_somatic+"PON-all_cns_varscan2.vcf.gz"
    output:
        output_somatic+"{sample}_somaticCall_varscan2_somatic_pon-filtered_strict_indel.vcf.gz"
    shell:
        """
        bcftools isec \
        -C \
        -c none \
        -w1 \
        {input.vcf} \
        {input.pon} | \
        bcftools filter \
        --include 'FORMAT/AD[0] = 0' \
        -o {output}
        """ 