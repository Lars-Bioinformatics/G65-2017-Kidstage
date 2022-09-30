__title__ = "Compute read counts in bam file for given positions"
__author__ = "Lars Andersen <larsmew@gmail.com>"
__date__ = "15/10/2022"
__version__ = "1.0"

import time, os, sys, glob

#########################################################
####                       Input                     ####
#########################################################
# Matched tumour-normal samples information
configfile: "/work/G65-2017-Kidstage/samples-matched.yaml"

# Explicit paths for external input files
ref = "/work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.fasta"
gnomad = "/work/sdukoldby/resources/hg38/af-only-gnomad.hg38.vcf.gz"
interval_list = "/work/sdukoldby/resources/hg38/MedExome_hg38_capture_targets.interval_list"
common_variants = "/work/sdukoldby/resources/hg38/small_exac_common_3.hg38.vcf.gz"

# Sample information
TUMORS = [tumor for tumor in config]
# print(PAIR)

map_qual = 30
base_qual = 30

# SAMPLES = "ECV2-29-blod_normal_tagseq-medexome.connor.recalibrated_for_pon.vcf"

#########################################################
####                   Input/Output                  ####
#########################################################
input_vcf = "mutect2_somatic_variants/vcf/vcf_filterFlag_PASS/"
output = "mutect2_somatic_variants/"


#########################################################
####                       Setup                     ####
#########################################################
# Timing
totim = time.time()
timeFormat = "%Y_%m_%d:%X" # year, month, day, time H:M:S

# Memory
mem = "-Xmx12g" # login nodes - be careful not running too many jobs at once!
# mem = "-Xmx24g" # slim nodes
# mem = "-Xmx32g" # Fat nodes
# mem = "-Xmx64g"

#########################################################
####  Define workflow start, stop and error actions  ####
#########################################################
# onstart:
#     shell("mkdir -p "+output_somatic)
    


#########################################################
####                  Run All Rules                  ####
#########################################################
'''
Rule all
'''
rule all:
    input:
        # [expand(output+"bam_readcount_output_q{map_qual}_b{base_qual}/{tumor}/combined_{tumor}_bam-readcount_VAF_q{map_qual}_b{base_qual}.txt", 
        #         tumor=tumor, map_qual=map_qual, base_qual=base_qual) for tumor in config],
        expand(output+"bam_readcount_output_q{map_qual}_b{base_qual}/bam-readcounts_vaf/{type}_bam-readcounts_vaf_boxplot_q{map_qual}_b{base_qual}.png", 
                tumor=TUMORS, map_qual=map_qual, base_qual=base_qual, type=["combined","snvs"]),
        # expand(output+"bam_readcount_output_q{map_qual}_b{base_qual}/{tumor}/{tumor}_bam-readcount_VAF_{type}_q{map_qual}_b{base_qual}.txt", 
            #    tumor="G65-T33B-4315_nimblegen-medexome_HYLKFDSXX", type="insertions", map_qual=map_qual, base_qual=base_qual),
        # [expand(output+"regions_bam_readcount/{tumor}_snvs_pass_fixed.bed", tumor=tumor) for tumor in config]
        # [expand(output+"bam_readcount_output_q20_b20/{tumor}/{type}/{sample}_{type}_readcounts_q20_b20.txt",
        #     tumor=tumor,
        #     type=type,
        #     # type="snvs",
        #     sample=config[tumor]["tumor"]) for tumor in config for type in ["snvs","insertions"]]
        # [expand(input_vcf+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag_PASSonly.vcf", 
        #                              normal=config[sample]["normal"],
        #                              tumor=sample) for sample in config]


# rule get_vcf_regions:
#     input:
#         vcf=input_vcf+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag_PASSonly.vcf.gz",
#     output:
#         bed=output+"regions_bam_readcount/{tumor}_vs_{normal}_pass_regions.bed"
#     shell:
#         """
#         zcat {input.vcf} | \
#         grep -v ^# | \
#         awk -F'\t' '{if (length($4)-1==0) printf($1"\t"$2"\t"$2+lenght($4)-1"\n")}' \
#         > {output.bed}
#         """

# rule decompress_vcf:
#     input:
#         vcf=input_vcf+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag_PASSonly.vcf.gz",
#     output:
#         vcf=input_vcf+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag_PASSonly.vcf",
#     shell:
#         """
#         bgzip -d -c {input.vcf} > {output.vcf}
#         """

rule vcf_norm_decompose_PASS_decompress:
    input:
        vcf=input_vcf+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag_PASSonly.vcf.gz",
    output:
        norm_vcf=input_vcf+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm_PASSonly.vcf",
    shell:
        """
        bcftools norm -m-any {input.vcf} | \
            bcftools norm -Ou --check-ref w -f {ref} | \
            bcftools filter -Ov -i 'FILTER="PASS"' > {output.norm_vcf}
        """

rule get_vcf_regions:
    input:
        vcf=lambda wildcards: expand(input_vcf+"{{tumor}}_vs_{normal}_somatic_mutect2_filterFlag_norm_PASSonly.vcf", 
                                     normal=config[wildcards.tumor]["normal"])
    output:
        bed_snvs=output+"regions_bam_readcount/bed/{tumor}_snvs_pass.bed",
        # bed_del=output+"regions_bam_readcount/bed/{tumor}_deletions_pass.bed",
        bed_ins=output+"regions_bam_readcount/bed/{tumor}_insertions_pass.bed"
    shell:
        # vcf2bed is from BEDOPS
        """
        vcf2bed -d --snvs < {input.vcf} > {output.bed_snvs}
        vcf2bed -d --insertions < {input.vcf} > {output.bed_ins}
        """
        # vcf2bed -d --deletions < {input.vcf} > {output.bed_del} ## Cannot be properly used at the moment

rule fix_vcf_regions:
    input:
        bed=output+"regions_bam_readcount/bed/{tumor}_{type}_pass.bed",
    output:
        bed=output+"regions_bam_readcount/{tumor}_{type}_pass_fixed.bed",
    shell:
        """
        awk -F'\t' '{{printf "%s\\t%s\\t%s\\t%s\\t%s\\n", $1, $2+1, $3, $6, $7}}' {input} > {output}
        """
        # awk -F'\t' '{printf "%s\t%s\t%s\t%s\t%s\n", $1, $2+1, $3, $6, $7}' {input} > {output}

rule bam_readcounts:
    input:
        bed=output+"regions_bam_readcount/{tumor}_snvs_pass_fixed.bed",
        bam="bam/{sample}.connor.recalibrated.bam"
    output:
        readcounts=output+"bam_readcount_output_q{map_qual}_b{base_qual}/{tumor}/snvs/{sample}_snvs_readcounts_q{map_qual}_b{base_qual}.txt"
    shell:
        """
        bam-readcount \
            --max-warnings 1 \
            -f {ref} \
            -q {map_qual} -b {base_qual} \
            -l {input.bed} \
            {input.bam} \
            > {output.readcounts}
        """
        # conda run --name bam_readcount \ - virker ikke for vcf2bed for whatever reason :/

rule bam_readcounts_insertions:
    input:
        bed=output+"regions_bam_readcount/{tumor}_insertions_pass_fixed.bed",
        bam="bam/{sample}.connor.recalibrated.bam"
    output:
        readcounts=output+"bam_readcount_output_q{map_qual}_b{base_qual}/{tumor}/insertions/{sample}_insertions_readcounts_q{map_qual}_b{base_qual}.txt"
    shell:
        """
        bam-readcount \
            --max-warnings 1 \
            --insertion-centric \
            -f {ref} \
            -q {map_qual} -b {base_qual} \
            -l {input.bed} \
            {input.bam} \
            > {output.readcounts}
        """

rule compute_readcounts_vaf:
    input:
        bed=output+"regions_bam_readcount/{tumor}_{type}_pass_fixed.bed",
        readcounts=lambda wildcards: expand(output+"bam_readcount_output_q{{map_qual}}_b{{base_qual}}/{{tumor}}/{{type}}/{sample}_{{type}}_readcounts_q{{map_qual}}_b{{base_qual}}.txt", 
            sample=config[wildcards.tumor]["tumor"])
    output:
        vaf=output+"bam_readcount_output_q{map_qual}_b{base_qual}/{tumor}/{tumor}_bam-readcount_VAF_{type}_q{map_qual}_b{base_qual}.txt"
    params:
        type = "{type}"
    script:
        "subscripts/bam-readcounts-VAF.py"

rule combine_vaf:
    input:
        vaf=expand(output+"bam_readcount_output_q{{map_qual}}_b{{base_qual}}/{{tumor}}/{{tumor}}_bam-readcount_VAF_{type}_q{{map_qual}}_b{{base_qual}}.txt", type=["snvs","insertions"])
    output:
        combined=output+"bam_readcount_output_q{map_qual}_b{base_qual}/{tumor}/combined_{tumor}_bam-readcount_VAF_q{map_qual}_b{base_qual}.txt"
    run:
        with open(output.combined,"w") as out:
            write_header = True
            for file in input:
                print(file)
                with open(file,"r") as f:
                    for linenum,line in enumerate(f):
                        # Write header only once
                        if linenum == 0:
                            if write_header == True:
                                write_header = False
                            else:
                                continue
                        # Write all other lines
                        out.write(line)

rule boxplot_and_vaf_tables_combined:
    input:
        vaf=expand(output+"bam_readcount_output_q{{map_qual}}_b{{base_qual}}/{tumor}/combined_{tumor}_bam-readcount_VAF_q{{map_qual}}_b{{base_qual}}.txt", tumor=TUMORS)
    output:
        table_long=output+"bam_readcount_output_q{map_qual}_b{base_qual}/bam-readcounts_vaf/combined_bam-readcounts_vaf_summarised_long_q{map_qual}_b{base_qual}.tsv",
        table_wide=output+"bam_readcount_output_q{map_qual}_b{base_qual}/bam-readcounts_vaf/combined_bam-readcounts_vaf_summarised_wide_q{map_qual}_b{base_qual}.tsv",
        boxplot=output+"bam_readcount_output_q{map_qual}_b{base_qual}/bam-readcounts_vaf/combined_bam-readcounts_vaf_boxplot_q{map_qual}_b{base_qual}.png"
    script:
        "subscripts/boxplot_and_vaf_tables.R"


rule boxplot_and_vaf_tables_snvs:
    input:
        vaf=expand(output+"bam_readcount_output_q{{map_qual}}_b{{base_qual}}/{tumor}/{tumor}_bam-readcount_VAF_snvs_q{{map_qual}}_b{{base_qual}}.txt", tumor=TUMORS)
    output:
        table_long=output+"bam_readcount_output_q{map_qual}_b{base_qual}/bam-readcounts_vaf/snvs_bam-readcounts_vaf_summarised_long_q{map_qual}_b{base_qual}.tsv",
        table_wide=output+"bam_readcount_output_q{map_qual}_b{base_qual}/bam-readcounts_vaf/snvs_bam-readcounts_vaf_summarised_wide_q{map_qual}_b{base_qual}.tsv",
        boxplot=output+"bam_readcount_output_q{map_qual}_b{base_qual}/bam-readcounts_vaf/snvs_bam-readcounts_vaf_boxplot_q{map_qual}_b{base_qual}.png"
    script:
        "subscripts/boxplot_and_vaf_tables.R"