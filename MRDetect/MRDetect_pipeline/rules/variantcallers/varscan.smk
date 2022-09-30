#########################################################
####                     VarScan                     ####
#########################################################

# Initial operations
varscan_output_dir = output_raw+"varscan/"
varscan_output_tmp = output_raw+"varscan/tmp/"
varscan_output_somatic = output_raw+"varscan/somatic_variants/"
onstart:
    shell("mkdir -p "+varscan_output_dir)
    shell("mkdir -p "+varscan_output_tmp)
    shell("mkdir -p "+varscan_output_somatic)

rule Varscan_Pileup:
    input:
        input_dir+"{sample}"+bamsuffix
    output:
        pileup=varscan_output_tmp+"{sample}_pileup.mpileup"
    resources:
        mem_mb=240000,
        cpus=threadsparallel
    threads:
        threadsparallel
    shell:
        """
        conda run --name samtools \
        samtools mpileup \
        --ignore-RG \
        --no-BAQ \
        --min-MQ 10 \
        --min-BQ 30 \
        --fasta-ref {ref} \
        --output {output.pileup} \
        {input}
        """

rule Varscan_SNPs_and_indels:
    input:
        normal=varscan_output_tmp+"{normal}_pileup.mpileup",
        tumor=varscan_output_tmp+"{tumor}_pileup.mpileup"
    output:
        snvs=varscan_output_somatic+"{tumor}_vs_{normal}_somatic_varscan_snvs.vcf",
        indel=varscan_output_somatic+"{tumor}_vs_{normal}_somatic_varscan_indels.vcf"
    resources:
        mem_mb=240000,
        cpus=1
    threads:
        1
    shell:
        """
        conda run --name varscan \
        varscan somatic \
        {input.normal} \
        {input.tumor} \
        --strand-filter 1 \
        --output-vcf 1 \
        --output-snp {output.snvs} \
        --output-indel {output.indel}
        """

checkpoint Varscan_copy:
    input:
        snvs=varscan_output_somatic+"{tumor}_vs_{normal}_somatic_varscan_snvs.vcf",
        indels=varscan_output_somatic+"{tumor}_vs_{normal}_somatic_varscan_indels.vcf"
    output:
        snvs=output_variants_notfiltered+"{tumor}_vs_{normal}/varscan_snvs.vcf",
        indels=output_variants_notfiltered+"{tumor}_vs_{normal}/varscan_indels.vcf"
    params:
        folder=output_variants_notfiltered+"{tumor}_vs_{normal}/"
    resources:
        mem_mb=12000,
        cpus=1
    threads:
        1
    shell:
        """
        mkdir -p {params.folder}
        cp {input.snvs} {output.snvs}
        cp {input.indels} {output.indels}
        """
