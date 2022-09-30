
#########################################################
####                       Setup                     ####
#########################################################
# onstart:
#     shell("mkdir -p "+output_somatic+"union_vcf/")


#########################################################
####                  Run All Rules                  ####
#########################################################
'''
Rule all
'''
# print(expand(output_somatic+"{sample}_somatic_mutect2_filterFlag.vcf.gz", sample=SAMPLES))

# rule all_pairs:
#     input:
#         # expand(output_somatic+"union_vcf/{sample}_mutect2_matched_filterFlag_norm_PASSonly_merged.vcf.gz", sample=SAMPLES)
#         # expand(output_somatic+"union_vcf/{sample}_mutect2_matchedAndJoint_filterFlag_norm_PASSonly_merged.vcf.gz", sample=SAMPLES)
#         expand(main_output+"union_vcf/{sample}_mutect2_matchedAndJoint_filterFlag_norm_PASSonly_merged.vcf.gz", sample=SAMPLES)


#1st pipe, splits multi-allelic calls into separate variant calls
#2nd pipe, left-aligns indels and issues warnings when the REF base in your VCF does not match the base in the supplied FASTA reference genome
rule vcf_joint_calling_norm:
    input:
        vcf="{output}vcf_filterFlag_PASS/{sample}_somatic_mutect2_filterFlag_PASSonly.vcf.gz"
    output:
        norm_vcf="{output}vcf_filterFlag_PASS_norm/{sample}_somatic_mutect2_filterFlag_norm_PASSonly.vcf.gz"
    shell:
        """
        bcftools norm -m-any {input.vcf} | \
            bcftools norm -Ou --check-ref w -f {ref} | \
            bcftools filter -Ou -i 'FILTER="PASS" | %FILTER="."' | \
            bcftools annotate -Oz -x 'ID' -I +'%CHROM:%POS:%REF:%ALT' > {output.norm_vcf}
        bcftools index -t {output.norm_vcf}
        """

# rule vcf_matched_norm:
#     input:
#         vcf=mutect2_matched_output+"vcf_filterFlag_PASS/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_PASSonly.vcf.gz"
#     output:
#         norm_vcf=mutect2_matched_output+"vcf_filterFlag_PASS_norm/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm_PASSonly.vcf.gz"
#     shell:
#         """
#         bcftools norm -m-any {input.vcf} | \
#             bcftools norm -Ou --check-ref w -f {ref} | \
#             bcftools filter -Ou -i 'FILTER="PASS"' | \
#             bcftools annotate -Oz -x 'ID' -I +'%CHROM:%POS:%REF:%ALT' > {output.norm_vcf}
#         bcftools index -t {output.norm_vcf}
#         """

# rule vcf_merge:
#     input:
#         norm_vcfs=lambda wildcards: expand(output_somatic+"union_vcf/norm_vcf/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm_PASSonly.vcf.gz", tumor=config[wildcards.sample]["tumor"], normal=config[wildcards.sample]["normal"])
#     output:
#         merged_vcf=output_somatic+"union_vcf/{sample}_mutect2_matched_filterFlag_norm_PASSonly_merged.vcf.gz"
#     shell:
#         """
#         bcftools merge -Oz -m none --force-samples {input.norm_vcfs} > {output.merged_vcf}
#         bcftools index -t {output.merged_vcf}
#         """


rule vcf_merge_matchedAndJoint:
    input:
        # norm_vcfs=lambda wildcards: expand(output_somatic+"union_vcf/norm_vcf/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm_PASSonly.vcf.gz", tumor=config[wildcards.sample]["tumor"], normal=config[wildcards.sample]["normal"]),
        norm_vcf_joint=mutect2_joint_calling_output+"vcf_filterFlag_PASS_norm/{sample}_somatic_mutect2_filterFlag_norm_PASSonly.vcf.gz",
        norm_vcfs_matched=lambda wildcards: expand(mutect2_matched_output+"vcf_filterFlag_PASS_norm/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm_PASSonly.vcf.gz",
               tumor=config[wildcards.sample]["tumor"],
               normal=config[wildcards.sample]["normal"])
    output:
        merged_vcf=main_output+"union_vcf/{sample}_mutect2_matchedAndJoint_filterFlag_norm_PASSonly_merged.vcf.gz"
    shell:
        """
        bcftools merge -Oz -m none --force-samples {input.norm_vcfs_matched} {input.norm_vcf_joint} > {output.merged_vcf}
        bcftools index -t {output.merged_vcf}
        """

# rule vcf_union:
#     input:
#         norm_vcfs=lambda wildcards: expand(output_somatic+"union_vcf/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm.vcf.gz", tumor=config[wildcards.sample]["tumor"], normal=config[wildcards.sample]["normal"])
#     output:
#         merged_vcf=output_somatic+"sites/sites.txt"
#     shell:
#         """
#         bcftools isec -Oz -n+1 -p sites {input.norm_vcfs}
#         """
