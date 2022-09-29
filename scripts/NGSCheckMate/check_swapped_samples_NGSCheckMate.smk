SAMPLES, = glob_wildcards("fastq/{sample}_R1.fastq.gz")

resource_path = "/work/sdukoldby/resources/hg38/"
ref = resource_path + "Homo_sapiens_assembly38.fasta"

# NGSCheckMate install dir
ncm_dir = "/work/G65-2017-Kidstage/NGSCheckMate"

rule all:
    input:
        "NGSCheckMate_output_bam/output_all.txt",
        # "NGSCheckMate_output/output_all.txt"


rule ngscheckmate_fastq:
    input:
        f1="fastq/{sample}_R1.fastq.gz",
        f2="fastq/{sample}_R2.fastq.gz",
        pt="NGSCheckMate/SNP/SNP.pt"
    output:
        vaf="NGSCheckMate_output/vaf/{sample}.vaf"
    shell:
        """
        ./NGSCheckMate/ngscheckmate_fastq -1 {input.f1} -2 {input.f2} {input.pt} > {output.vaf}
        """

rule vaf_ncm:
    input:
        vaf=expand("NGSCheckMate_output/vaf/{sample}.vaf", sample=SAMPLES)
    output:
        "NGSCheckMate_output/output_all.txt"
    params:
        in_dir=directory("NGSCheckMate_output/vaf/"),
        out_dir=directory("NGSCheckMate_output"),
        prefix="output"
    shell:
        """
        python2 {ncm_dir}/vaf_ncm.py -f -I {params.in_dir} -O {params.out_dir} {params.prefix}
        """

rule vcf_from_bam:
    input:
        bam="Connor/bam/{sample}.connor.recalibrated.bam",
        bed=ncm_dir+"/SNP/SNP_GRCh38_hg38_wChr.bed",
    output:
        vcf="NGSCheckMate_output_bam/vcf/{sample}.vcf"
    shell:
        """
        samtools mpileup -I -uf {ref} -l {input.bed} {input.bam} | bcftools call -c - > {output.vcf}
        """

rule ngscheckmate_bam:
    input:
        vcf=expand("NGSCheckMate_output_bam/vcf/{sample}.vcf", sample=SAMPLES),
        bed=ncm_dir+"/SNP/SNP_GRCh38_hg38_wChr.bed",
        # sampleList="NGSCheckMate_output_bam/NGSCheckMate_vcf_sampleList.txt"
    output:
        "NGSCheckMate_output_bam/output_all.txt"
    params:
        in_dir=directory("NGSCheckMate_output_bam/vcf"),
        out_dir=directory("NGSCheckMate_output_bam")
    shell:
        """
        python2 {ncm_dir}/ncm.py -V -d {params.in_dir} -bed {input.bed} -O {params.out_dir}
        """
