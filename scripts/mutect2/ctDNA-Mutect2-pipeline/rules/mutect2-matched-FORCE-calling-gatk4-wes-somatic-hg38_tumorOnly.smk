# onstart:
#     # shell("mkdir -p "+"{output}split/{{vcf,bam,f1r2}}_split")
#     shell("mkdir -p "+"{output}vcf_{{files,filterFlag,filterFlag_PASS,filterFlag_norm_PASS}}")
#     # shell("mkdir -p "+"{output}stats/")
#     shell("mkdir -p "+"{output}realigned_bam/")


#########################################################
####                  Run All Rules                  ####
#########################################################
'''
Rule all
'''
# print(expand("{output}{sample}_somatic_mutect2_filterFlag.vcf.gz", sample=SAMPLES))
# print(expand("{output}vcf_filterFlag_norm_decomposed/{sample}_somatic_mutect2_filterFlag_norm.vcf.gz", sample=SAMPLES))
# sys.exit()

# rule all_pairs:
#     input:
#         # expand("{output}vcf_filterFlag_PASS/{sample}_somatic_mutect2_filterFlag_PASSonly.vcf.gz", sample=SAMPLES)
#         expand("{output}vcf_final/{sample}_somatic_mutect2_filterFlag_norm_final.vcf.gz", sample=SAMPLES)
        # expand("{output}vcf_filterFlag_norm/{sample}_somatic_mutect2_filterFlag_norm.vcf.gz", sample=SAMPLES)
        # expand("{output}vcf_filterFlag_norm_PASS/{sample}_somatic_mutect2_filterFlag_norm_PASSonly.vcf.gz", sample=SAMPLES)
        # expand("{output}vcf_filterFlag_norm_decomposed/{sample}_somatic_mutect2_filterFlag_norm.vcf.gz", sample=SAMPLES)
        # expand("{output}vcf_files/{sample}_somatic_mutect2.vcf.gz", sample=SAMPLES),
        # expand("{output}{sample}_somatic_mutect2_filterFlag.vcf.gz", sample=SAMPLES),
        # expand("{output}{sample}_somatic_mutect2.bam", sample=SAMPLES)
        # expand("{output}{sample}_somatic_mutect2.vcf", sample=SAMPLES),
        # expand("{output}{sample}_merged_contamination.table", sample=SAMPLES)


##########################################################
####  Call Somatic Variants using Mutect2 on matched  ####
####   Tumor-Normal samples on per chromesome basis   ####
##########################################################
# Assumes PON already calculated in other mutect2 script.
'''
Mutect2 on matched Tumor-Normal samples
'''
rule Mutect2_matched_calling_forced_tumorOnly:
    input:
        # normal="bam/{normal}.connor.recalibrated.bam",
        tumor="bam/{tumor}.connor.recalibrated.bam",
        # pon=pon_location+"pon.vcf.gz",
        # intervals=target_regions+"{interval}.interval_list"
        force_intervals=lambda wildcards: expand(main_output+"union_vcf/{patient}"+force_intervals_suffix, patient=wildcards.tumor[:5]+"_thruplex-tag-seq-hv-nimblegen-medexome_HVFFJDSXY")
    output:
        vcf="{output}vcf_files_tumorOnly/{tumor}_somatic_mutect2-forced.vcf.gz",
        idx="{output}vcf_files_tumorOnly/{tumor}_somatic_mutect2-forced.vcf.gz.tbi",
        # realigned_bam="{output}realigned_bam/{tumor}_somatic_mutect2-forced.bam",
        vcf_stats="{output}vcf_files_tumorOnly/{tumor}_somatic_mutect2-forced.vcf.gz.stats",
        f1r2="{output}f1r2_tumorOnly/{tumor}_f1r2.tar.gz"
    resources: cpus=4, mem_mb=24000
    shell:
        """
        gatk --java-options -Xmx{resources.mem_mb}g Mutect2 \
        -R {ref} \
        -I {input.tumor} \
        -pon {pon_1000g} \
        --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --af-of-alleles-not-in-resource 0.0000025 \
        --germline-resource {gnomad} \
        --native-pair-hmm-threads {resources.cpus} \
        --f1r2-tar-gz {output.f1r2} \
        -L {input.force_intervals} \
        -alleles {input.force_intervals} \
        -O {output.vcf}
        """
        # -bamout {output.realigned_bam} \
        # -I {input.normal} \
        # -normal {wildcards.normal} \

rule FilterMutectCalls_forced_tumorOnly:
    input:
        vcf="{output}vcf_files_tumorOnly/{sample}_somatic_mutect2-forced.vcf.gz",
        idx="{output}vcf_files_tumorOnly/{sample}_somatic_mutect2-forced.vcf.gz.tbi",
        # bam="{output}realigned_bam/{sample}_somatic_mutect2-forced.bam",
        stats="{output}vcf_files_tumorOnly/{sample}_somatic_mutect2-forced.vcf.gz.stats",
        # contamination=mutect2_matched_output+"contaminations_tables/{sample}_merged_contamination.table",
        # read_orientation="{output}read_orientation_models/{sample}_read-orientation-model.tar.gz"
    output:
        vcf="{output}vcf_filterFlag_tumorOnly/{sample}_somatic_mutect2-forced_filterFlag.vcf.gz",
        tbi="{output}vcf_filterFlag_tumorOnly/{sample}_somatic_mutect2-forced_filterFlag.vcf.gz.tbi",
    shell:
        """
        gatk --java-options {mem} FilterMutectCalls \
        -R {ref} \
        -V {input.vcf} \
        --stats {input.stats} \
        -L {interval_list} \
        -O {output.vcf}
        """
        # --contamination-table {input.contamination} \
        # --orientation-bias-artifact-priors {input.read_orientation} \

rule vcf_norm_decompose_bcftools_tumorOnly:
    input:
        vcf="{output}vcf_filterFlag_tumorOnly/{sample}_somatic_mutect2-forced_filterFlag.vcf.gz"
    output:
        norm_vcf="{output}vcf_filterFlag_norm_tumorOnly/{sample}_somatic_mutect2-forced_filterFlag_norm.vcf.gz"
    shell:
        """
        bcftools norm -m-any {input.vcf} | \
            bcftools norm -Ou --check-ref w -f {ref} | \
            bcftools annotate -Oz -x 'ID' -I +'%CHROM:%POS:%REF:%ALT' > {output.norm_vcf}
        bcftools index -t {output.norm_vcf}
        """

rule intersect_forcedVCF_mergedVCF_tumorOnly:
    input:
        norm_vcf="{output}vcf_filterFlag_norm_tumorOnly/{tumor}_somatic_mutect2-forced_filterFlag_norm.vcf.gz",
        force_intervals=lambda wildcards: expand(main_output+"union_vcf/{patient}"+force_intervals_suffix, patient=wildcards.tumor[:5]+"_thruplex-tag-seq-hv-nimblegen-medexome_HVFFJDSXY")
    output:
        final_vcf="{output}vcf_final_tumorOnly/{tumor}__somatic_mutect2-forced_filterFlag_norm_final.vcf.gz"
    shell:
        """
        bcftools isec -n+2 -w1 -Oz {input.norm_vcf} {input.force_intervals} > {output.final_vcf}
        bcftools index -t {output.final_vcf}
        """

# rule vcf_norm_decompose:
#     input:
#         vcf="{output}vcf_filterFlag/{sample}_somatic_mutect2_filterFlag.vcf.gz"
#     output:
#         norm_vcf="{output}vcf_filterFlag_norm_decomposed/{sample}_somatic_mutect2_filterFlag_norm.vcf.gz"
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
#         norm_vcf="{output}vcf_filterFlag_norm_decomposed/{sample}_somatic_mutect2_filterFlag_norm.vcf.gz"
#     output:
#         vcf_pass="{output}vcf_filterFlag_norm_decomposed/{sample}_somatic_mutect2_filterFlag_norm_PASSonly.vcf.gz"
#     shell:
#         """
#         zcat {input.norm_vcf} | egrep '^#|PASS' | bgzip -c > {output.vcf_pass}
#         tabix -p vcf {output.vcf_pass}
#         """