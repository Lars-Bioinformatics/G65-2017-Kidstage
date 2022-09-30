#########################################################
####                 Somatic Sniper                  ####
#########################################################

# Initial operations
somaticsniper_output_dir = output_raw+"somaticsniper/"
somaticsniper_output_tmp = output_raw+"somaticsniper/tmp/"
somaticsniper_output_somatic = output_raw+"somaticsniper/somatic_variants/"
onstart:
    shell("mkdir -p "+somaticsniper_output_dir)
    shell("mkdir -p "+somaticsniper_output_tmp)
    shell("mkdir -p "+somaticsniper_output_somatic)

rule SomaticSniper_call:
    input:
        normal=input_dir+"{normal}"+bamsuffix,
        tumor=input_dir+"{tumor}"+bamsuffix
    output:
        vcf=somaticsniper_output_tmp+"{tumor}_vs_{normal}_somatic_somaticsniper_snvs.vcf"
    resources:
        mem_mb=240000,
        cpus=1
    threads:
        1
    shell:
        """
        conda run --name somatic-sniper \
        bam-somaticsniper \
        -Q 40 \
        -G \
        -L \
        -f {ref} \
        -F vcf \
        {input.tumor} \
        {input.normal} \
        {output.vcf}
        """

rule SomaticSniper_pileup:
    input:
        tumor=input_dir+"{tumor}"+bamsuffix
    output:
        pileup=somaticsniper_output_tmp+"{tumor}.pileup"
    resources:
        mem_mb=150000,
        cpus=1
    threads: 
        1
    shell:
        """
        conda run --name somatic-sniper \
        samtools pileup \
        -vci \
        -f {ref} \
        {input.tumor} > {output.pileup}
        """

rule SomaticSniper_pileup_filter:
    input:
        pileup=somaticsniper_output_tmp+"{tumor}.pileup"
    output:
        pileup=somaticsniper_output_tmp+"{tumor}_filtered.pileup"
    resources:
        mem_mb=12000,
        cpus=1
    threads:
        1
    shell:
        """
        conda run scripts/samtools/varfilter.py \
        -Q 40 {input.pileup} > {output.pileup}
        """

rule SomaticSniper_snpfilter:
    input:
        vcf=somaticsniper_output_tmp+"{tumor}_vs_{normal}_somatic_somaticsniper_snvs.vcf",
        pileup=somaticsniper_output_tmp+"{tumor}_filtered.pileup"
    output:
        snpfilter=somaticsniper_output_tmp+"{tumor}_vs_{normal}_somatic_somaticsniper_snvs.vcf.SNPfilter"
    resources:
        mem_mb=12000,
        cpus=1
    threads:
        1
    shell:
        """
        perl scripts/SomaticSniper/snpfilter.pl \
        --snp-file {input.vcf} \
        --indel-file {input.pileup}
        """

rule SomaticSniper_prepare_for_readcount:
    input:
        snpfilter=somaticsniper_output_tmp+"{tumor}_vs_{normal}_somatic_somaticsniper_snvs.vcf.SNPfilter"
    output:
        snpfilter=somaticsniper_output_tmp+"{tumor}_vs_{normal}_somatic_somaticsniper_snvs.vcf.SNPfilter.pos"
    resources:
        mem_mb=12000,
        cpus=1
    threads:
        1
    shell:
        """
        perl scripts/SomaticSniper/prepare_for_readcount.pl \
        --snp-file {input.snpfilter}
        """

rule SomaticSniper_readcount:
    input:
        snpfilter=somaticsniper_output_tmp+"{tumor}_vs_{normal}_somatic_somaticsniper_snvs.vcf.SNPfilter.pos",
        tumor=input_dir+"{tumor}"+bamsuffix
    output:
        rc=somaticsniper_output_tmp+"{tumor}_vs_{normal}_somatic_somaticsniper_snvs.vcf.rc"
    resources:
        mem_mb=150000,
        cpus=1
    threads:
        1
    shell:
        """
        conda run --name readcount \
        bam-readcount \
        -b 40 \
        -f {ref} \
        -l {input.snpfilter} \
        {input.tumor} > {output.rc}
        """

rule SomaticSniper_fpfilter:
    input:
        snpfilter=somaticsniper_output_tmp+"{tumor}_vs_{normal}_somatic_somaticsniper_snvs.vcf.SNPfilter",
        rc=somaticsniper_output_tmp+"{tumor}_vs_{normal}_somatic_somaticsniper_snvs.vcf.rc"
    output:
        snpfilter=somaticsniper_output_tmp+"{tumor}_vs_{normal}_somatic_somaticsniper_snvs.vcf.SNPfilter.fp_pass"
    resources:
        mem_mb=100000,
        cpus=1
    threads:
        1
    shell:
        """
        perl scripts/SomaticSniper/fpfilter.pl \
        --snp-file {input.snpfilter} \
        --readcount-file {input.rc}
        """

rule SomaticSniper_highconfidence:
    input:
        snpfilter=somaticsniper_output_tmp+"{tumor}_vs_{normal}_somatic_somaticsniper_snvs.vcf.SNPfilter.fp_pass"
    output:
        snpfilter=somaticsniper_output_tmp+"{tumor}_vs_{normal}_somatic_somaticsniper_snvs.vcf.SNPfilter.fp_pass.hc"
    resources:
        mem_mb=12000,
        cpus=1
    threads:
        1
    shell:
        """
        perl scripts/SomaticSniper/highconfidence.pl \
        --snp-file {input.snpfilter}
        """

checkpoint SomaticSniper_copy:
    input:
        snvs=somaticsniper_output_tmp+"{tumor}_vs_{normal}_somatic_somaticsniper_snvs.vcf.SNPfilter.fp_pass"
    output:
        snvs=output_variants_notfiltered+"{tumor}_vs_{normal}/somaticsniper_snvs.vcf"
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
        awk '
        {{OFS = "\t"}} 
        /^#/ {{print}}
        !/^#/ {{print $1,$2,$3,$4,$5,$6,"PASS",$8,$9,$10,$11}}
        ' {input.snvs} > {output.snvs}
        """