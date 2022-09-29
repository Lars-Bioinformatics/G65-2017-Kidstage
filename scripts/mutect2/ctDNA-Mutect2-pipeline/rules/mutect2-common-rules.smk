rule merge_somatic_vcf:
    input:
        vcf_subfile=expand("{{output}}split/vcf_split/{{sample}}_somatic_mutect2__{interval}__split.vcf.gz", interval=INTERVALS)
    output:
        vcf="{output}vcf_files/{sample}_somatic_mutect2.vcf.gz",
        idx="{output}vcf_files/{sample}_somatic_mutect2.vcf.gz.tbi",
    params:
        vcf_subfile=expand("-I "+"{{output}}split/vcf_split/{{sample}}_somatic_mutect2__{interval}__split.vcf.gz", interval=INTERVALS)
    shell:
        """
        gatk --java-options {mem} GatherVcfs \
        {params.vcf_subfile} \
        -O {output.vcf}

        tabix -p vcf {output.vcf}
        """

## FOR DEBUGGING
rule merge_mutect2_bam:
  input:
        bam_subfile=expand("{{output}}split/bam_split/{{sample}}__{interval}_mutect2.bam", interval=INTERVALS)
  output:
        bam="{output}realigned_bam/{sample}_somatic_mutect2.bam"
  params:
        bam_subfile=expand("-I "+"{{output}}split/bam_split/{{sample}}__{interval}_mutect2.bam", interval=INTERVALS)
  shell:
        """
        gatk --java-options {mem} MergeSamFiles \
        {params.bam_subfile} \
        --VALIDATION_STRINGENCY LENIENT \
        --USE_THREADING true \
        -O {output.bam} \
        --CREATE_INDEX true
        """

rule merge_somatic_vcf_stats:
    input:
        vcf_subfile=expand("{{output}}split/vcf_split/{{sample}}_somatic_mutect2__{interval}__split.vcf.gz.stats", interval=INTERVALS)
    output:
        vcf_stats="{output}vcf_files/{sample}_somatic_mutect2.vcf.gz.stats",
    params:
        vcf_subfile=expand("-stats "+"{{output}}split/vcf_split/{{sample}}_somatic_mutect2__{interval}__split.vcf.gz.stats", interval=INTERVALS)
    shell:
        """
        gatk --java-options {mem} MergeMutectStats \
        {params.vcf_subfile} \
        -O {output.vcf_stats}
        """

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
        f1r2=expand("{{output}}split/f1r2_split/{{sample}}_f1r2__{interval}__split.tar.gz", interval=INTERVALS)
    output:
        f1r2_model="{output}read_orientation_models/{sample}_read-orientation-model.tar.gz"
    params:
        f1r2=expand("-I "+"{{output}}split/f1r2_split/{{sample}}_f1r2__{interval}__split.tar.gz", interval=INTERVALS)
    shell:
        """
        gatk LearnReadOrientationModel \
        {params.f1r2} \
        -O {output}
        """
        # -I {input} \

#########################################################
####           Create Contamination table            ####
#########################################################
### Sæt sammen til en metode, som tager højde for med og uden normal... (Se gatk4 scripts)
rule GetPileupSummaries:
    input:
        bam="bam/{sample}.connor.recalibrated.bam"
    output:
        pileup="{output}pileups/{sample}_pileup.table"
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
        normal="{output}pileups/{normal}_pileup.table",
        tumor="{output}pileups/{tumor}_pileup.table"
    output:
        contamination="{output}contaminations_tables/{tumor}_vs_{normal}_contamination.table"
    shell:
        """
        gatk --java-options {mem} CalculateContamination \
        -I {input.tumor} \
        -matched {input.normal} \
        -O {output}
        """

rule JoinContaminationTables:
    input:
        cont_tables=lambda wildcards: expand("{{output}}contaminations_tables/{tumor}_vs_{normal}_contamination.table", normal=config[wildcards.sample]["normal"], tumor=config[wildcards.sample]["tumor"])
    output:
        joined_table="{output}contaminations_tables/{sample}_merged_contamination.table"
    shell:
        """
        echo -e 'sample\tcontamination\terror' > {output};
        awk -F'\t' 'FNR == 2' {input} >> {output}
        """
    
rule JoinContaminationTables_matched:
    input:
        cont_tables="{output}/contaminations_tables/{tumor}_vs_{normal}_contamination.table"
    output:
        joined_table="{output}/contaminations_tables/{tumor}_vs_{normal}_merged_contamination.table"
    shell:
        """
        echo -e 'sample\tcontamination\terror' > {output};
        awk -F'\t' 'FNR == 2' {input} >> {output}
        """


#########################################################
####              Filter Mutect2 Calls               ####
#########################################################
rule FilterMutectCalls:
    input:
        vcf="{output}vcf_files/{sample}_somatic_mutect2.vcf.gz",
        idx="{output}vcf_files/{sample}_somatic_mutect2.vcf.gz.tbi",
        bam="{output}realigned_bam/{sample}_somatic_mutect2.bam",
        stats="{output}vcf_files/{sample}_somatic_mutect2.vcf.gz.stats",
        contamination="{output}contaminations_tables/{sample}_merged_contamination.table",
        read_orientation="{output}read_orientation_models/{sample}_read-orientation-model.tar.gz"
    output:
        vcf="{output}vcf_filterFlag/{sample}_somatic_mutect2_filterFlag.vcf.gz",
        tbi="{output}vcf_filterFlag/{sample}_somatic_mutect2_filterFlag.vcf.gz.tbi",
    shell:
        """
        gatk --java-options {mem} FilterMutectCalls \
        -R {ref} \
        -V {input.vcf} \
        --contamination-table {input.contamination} \
        --orientation-bias-artifact-priors {input.read_orientation} \
        --stats {input.stats} \
        -L {interval_list} \
        -O {output.vcf}
        """

rule vcf_pass_only:
    input:
        vcf="{output}vcf_filterFlag/{sample}_somatic_mutect2_filterFlag.vcf.gz",
    output:
        vcf_pass="{output}vcf_filterFlag_PASS/{sample}_somatic_mutect2_filterFlag_PASSonly.vcf.gz",
    shell:
        """
        bcftools view -i "%FILTER='PASS' | %FILTER='.'" {input.vcf} | bgzip -c > {output.vcf_pass}
        tabix -p vcf {output.vcf_pass}
        """

rule vcf_norm_decompose:
    input:
        vcf="{output}vcf_filterFlag/{sample}_somatic_mutect2_filterFlag.vcf.gz"
    output:
        norm_vcf="{output}vcf_filterFlag_norm_decomposed/{sample}_somatic_mutect2_filterFlag_norm.vcf.gz"
    shell:
        """
        zless {input.vcf} \
        | vt decompose -s - \
        | vt normalize -r {ref} - \
        | bgzip -c > {output.norm_vcf}
        """
        # vcf index: tabix -p vcf {output.norm_vcf}
        # Alternative to vt
        # bcftools norm -m-any {input.vcf} > {output.vcf}

rule vcf_norm_pass_only:
    input:
        norm_vcf="{output}vcf_filterFlag_norm_decomposed/{sample}_somatic_mutect2_filterFlag_norm.vcf.gz"
    output:
        vcf_pass="{output}vcf_filterFlag_norm_decomposed/{sample}_somatic_mutect2_filterFlag_norm_PASSonly.vcf.gz"
    shell:
        """
        bcftools view -i "%FILTER='PASS' | %FILTER='.'" {input.norm_vcf} | bgzip -c > {output.vcf_pass}
        tabix -p vcf {output.vcf_pass}
        """