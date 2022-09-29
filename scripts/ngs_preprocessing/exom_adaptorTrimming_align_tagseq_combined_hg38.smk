# Run script from folder with fastq files.

from os import getcwd
import time

# SAMPLES, = glob_wildcards("{sample}_R1_001.fastq.gz")
SAMPLES, = glob_wildcards("fastq/{sample}_R1.fastq.gz")
#SAMPLES = "PC1-10-blod_normal_tagseq-medexome"
#SAMPLES = "PC1-14-recidiv_tumor_tagseq-medexome"
FAMNAME = getcwd().rsplit("/",1)[1]

resource_path = "/work/sdukoldby/resources/hg38/"
ref = resource_path + "Homo_sapiens_assembly38.fasta"
dbsnp= resource_path + "Homo_sapiens_assembly38.dbsnp138.vcf"
mills_1000g=resource_path + "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
phase1_1000g=resource_path + "1000G_phase1.snps.high_confidence.hg38.vcf.gz"
bed = resource_path + "MedExome_hg38_capture_targets.bed"
interval_list = resource_path + "MedExome_hg38_capture_targets.interval_list"

### Can be run without bed file

mem = 50
timeFormat = time.strftime("%Y_%m_%d-%X")

rule all:
  input:
    # expand("MarkDuplicates/bam/{sample}.recalibrated.bam", sample=SAMPLES),
    # expand("bam/{sample}.connor.recalibrated.bam", sample=SAMPLES),
    "logs/exom_align_tagseq_combined_hg38-version-log_"+timeFormat+".txt"


rule fastqc:
    input:
        "fastq/{sample}_R1.fastq.gz",
        "fastq/{sample}_R2.fastq.gz",
    output:
        html_1 = "quality_control/fastqc/{sample}_R1_fastqc.html",
        zip_1 = "quality_control/fastqc/{sample}_R1_fastqc.zip",
        html_2 = "quality_control/fastqc/{sample}_R2_fastqc.html",
        zip_2 = "quality_control/fastqc/{sample}_R2_fastqc.zip",
    shell:
        """
        fastqc {input} -o "quality_control/fastqc/"
        """

# Check for correct adaptors:
# zcat f1.fastq.gz | head -4000000 | grep AGATCGGAAGAGCACACGTCTGAACTCCAGTCA | wc -l
# zcat f2.fastq.gz | head -4000000 | grep AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT | wc -l
rule adapter_trimming:
    input:
        f1="fastq/{sampleid}_{protocol}_{flowcell}_R1.fastq.gz",
        f2="fastq/{sampleid}_{protocol}_{flowcell}_R2.fastq.gz",
    output:
        f1="fastq_trimmed/{sampleid}_{protocol}_{flowcell}-trimmed-pair1.fastq.gz",
        f2="fastq_trimmed/{sampleid}_{protocol}_{flowcell}-trimmed-pair2.fastq.gz"
    params:
        fastq_name="fastq_trimmed/{sampleid}_{protocol}_{flowcell}"
    resources:
        cpus=7, mem_mb=42000
    shell:
        """
        skewer \
        -x AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -y AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -t {resources.cpus} \
        --compress \
        -o {params.fastq_name} \
        {input}
        """
rule bwa_mem:
    input:
        "fastq_trimmed/{sampleid}_{protocol}_{flowcell}-trimmed-pair1.fastq.gz",
        "fastq_trimmed/{sampleid}_{protocol}_{flowcell}-trimmed-pair2.fastq.gz"
    output:
        bam = temp("{sampleid}_{protocol}_{flowcell}.sorted.bam"),
        bai = temp("{sampleid}_{protocol}_{flowcell}.sorted.bai")
    resources: cpus=21, mem_mb=126000
    log:
        "logs/bwa_logs/{sampleid}_{protocol}_{flowcell}.bwa_mem.log"
    params:
        rgid = "{sampleid}_{protocol}_{flowcell}",
        rglb = "{protocol}",
        rgsm = "{sampleid}_{protocol}_{flowcell}",
        rgpl = "illumina",
        rgpu = "{flowcell}"
    shell:
        """
        bwa mem {ref} {input} \
        -R \"@RG\\tID:{params.rgid}\\tLB:{params.rglb}\\tSM:{params.rgsm}\\tPL:{params.rgpl}\\tPU:{params.rgpu}\" \
        -M \
        -t {resources.cpus} | \
        gatk --java-options -Xmx{mem}G SortSam \
        --INPUT /dev/stdin \
        --OUTPUT {output.bam} \
        --VALIDATION_STRINGENCY LENIENT \
        --SORT_ORDER coordinate \
        --TMP_DIR tmp \
        --CREATE_INDEX TRUE
        """

#####################################
####       Markdup pipeline      ####
#####################################
rule markdup:
    input:
        bam = "{sample}.sorted.bam",
        bai = "{sample}.sorted.bai"
    output:
        bam = temp("MarkDuplicates/bam/{sample}.marked.bam"),
        bai = temp("MarkDuplicates/bam/{sample}.marked.bai"),
        metrics = "MarkDuplicates/quality_control/{sample}_markdup_metrics.txt"
    resources: cpus=2, mem_mb=12000
    log:
        "logs/markdup_logs/{sample}.markdup.log"
    shell:
        """
        gatk --java-options -Xmx{mem}G MarkDuplicates \
        --INPUT {input.bam} \
        --OUTPUT {output.bam} \
        --METRICS_FILE {output.metrics} \
        --REMOVE_DUPLICATES FALSE \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
        --CREATE_INDEX TRUE
        """

# OPTICAL_DUPLICATE_PIXEL_DISTANCE:
# Should be set to 2500 for patterned flowcells like Illumina Hiseq 3000/4000 and Novaseq 6000
# and set to 100 for unpatterned flowcells like HiSeq 2500

rule BaseRecalibrator_markdup:
    input:
        bam = "MarkDuplicates/bam/{sample}.marked.bam",
        bai = "MarkDuplicates/bam/{sample}.marked.bai",
    output:
        "MarkDuplicates/quality_control/{sample}_pre_recalibration.grp"
    shell:
        """
        gatk --java-options -Xmx{mem}G BaseRecalibrator \
        -R {ref} \
        -I {input.bam} \
        --known-sites {dbsnp} \
        --known-sites {mills_1000g} \
        --known-sites {phase1_1000g} \
        -O {output}
        """

rule ApplyBQSR_markdup:
    input:
        bam = "MarkDuplicates/bam/{sample}.marked.bam",
        bai = "MarkDuplicates/bam/{sample}.marked.bai",
        recal = "MarkDuplicates/quality_control/{sample}_pre_recalibration.grp",
    output:
        bam = "MarkDuplicates/bam/{sample}.recalibrated.bam",
        bai = "MarkDuplicates/bam/{sample}.recalibrated.bai"
    shell:
        """
        gatk --java-options -Xmx{mem}G ApplyBQSR \
        -R {ref} \
        -I {input.bam} \
        --bqsr {input.recal} \
        -O {output.bam}
        """
#-L {input.bed}
#bed = {bed}

rule post_recalibration_table_markdup:
    input:
        bam = "MarkDuplicates/bam/{sample}.recalibrated.bam",
        bai = "MarkDuplicates/bam/{sample}.recalibrated.bai",
    output:
        "MarkDuplicates/quality_control/{sample}_post_recalibration.grp"
    shell:
        """
        gatk --java-options -Xmx{mem}G BaseRecalibrator \
        -R {ref} \
        -I {input.bam} \
        --known-sites {dbsnp} \
        --known-sites {mills_1000g} \
        --known-sites {phase1_1000g} \
        -O {output}
        """

#-L {input.bed}
#bed = {bed}

rule CollectHsMetrics_markdup:
    input:
        bam = "MarkDuplicates/bam/{sample}.recalibrated.bam",
        bai = "MarkDuplicates/bam/{sample}.recalibrated.bai",
        interval = {interval_list}
    output:
        "MarkDuplicates/quality_control/{sample}_capture-targets_HsMetrics.txt"
    shell:
        """
        gatk --java-options -Xmx{mem}G CollectHsMetrics \
        --INPUT {input.bam} \
        --REFERENCE_SEQUENCE {ref} \
        --OUTPUT {output} \
        --BAIT_INTERVALS {input.interval} \
        --TARGET_INTERVALS 	{input.interval}
        """

rule CollectInsertSizeMetrics_markdup:
    input:
        bam = "MarkDuplicates/bam/{sample}.recalibrated.bam"
    output:
        hist = "MarkDuplicates/quality_control/{sample}.insert_size_histogram.pdf",
        metrics = "MarkDuplicates/quality_control/{sample}.insert_size_metrics.txt"
    shell:
        """
        gatk --java-options -Xmx{mem}G CollectInsertSizeMetrics \
        -I {input.bam} \
        -H {output.hist} \
        -O {output.metrics}
        """


rule qualimap_markdup:
    input:
        "MarkDuplicates/bam/{sample}.recalibrated.bam"
    output:
        html = "MarkDuplicates/quality_control/{sample}/qualimapReport.html"
    params:
        outdir = directory("MarkDuplicates/quality_control/{sample}")
    resources: cpus=4, mem=30000
    shell:
        """
        qualimap bamqc \
        -bam {input} \
        -nt 6 \
        -c \
        -sd \
        -gff {bed} \
        -outdir {params.outdir} \
        --java-mem-size=20G
        """



rule multiqc_markdup:
    input:
        expand("quality_control/fastqc/{sample}_R1_fastqc.html", sample=SAMPLES),
        expand("quality_control/fastqc/{sample}_R1_fastqc.zip", sample=SAMPLES),
        expand("quality_control/fastqc/{sample}_R2_fastqc.html", sample=SAMPLES),
        expand("quality_control/fastqc/{sample}_R2_fastqc.zip", sample=SAMPLES),
        expand("MarkDuplicates/quality_control/{sample}_pre_recalibration.grp", sample=SAMPLES),
        expand("MarkDuplicates/quality_control/{sample}/qualimapReport.html", sample=SAMPLES),
        expand("MarkDuplicates/quality_control/{sample}_capture-targets_HsMetrics.txt", sample=SAMPLES),
        expand("MarkDuplicates/quality_control/{sample}_post_recalibration.grp", sample=SAMPLES),
        expand("MarkDuplicates/quality_control/{sample}_markdup_metrics.txt", sample=SAMPLES)
    output:
        html = "MarkDuplicates/quality_control/{famname}_multiqc_report.html"
    params:
        indir1 = directory("quality_control/fastqc/"),
        indir2 = directory("MarkDuplicates/quality_control/"),
        config = resource_path + "multiqc_config.yaml"
    shell:
        """
        multiqc {params.indir1} {params.indir2} \
        -n {output.html} \
        -c {params.config}
        """

#####################################
####       Connor pipeline       ####
#####################################

rule connor:
    input:
        bam = "{sample}.sorted.bam",
        bai = "{sample}.sorted.bai"
    output:
        bam = temp("Connor/bam/{sample}.connor.dedup.bam"),
        bai = temp("Connor/bam/{sample}.connor.dedup.bam.bai")
    log:
        "Connor/quality_control/{sample}.connor_dedup.log"
    shell:
         """
         connor {input.bam} {output.bam} \
         -s 0 \
         --force \
         --log_file={log}
         """

#--force er nødvendig for at få connor til at køre på mergede fastq-filer hvor data har forskellige read lengths.


rule BaseRecalibrator_connor:
    input:
        bam = "Connor/bam/{sample}.connor.dedup.bam",
        bai = "Connor/bam/{sample}.connor.dedup.bam.bai",
    output:
        "Connor/quality_control/{sample}_pre_recalibration.grp"
    shell:
        """
        gatk --java-options -Xmx{mem}G BaseRecalibrator \
        -R {ref} \
        -I {input.bam} \
        --known-sites {dbsnp} \
        --known-sites {mills_1000g} \
        --known-sites {phase1_1000g} \
        -O {output}
        """
#        -L {input.bed}
#bed = {bed}
### Can be run without bed file


rule ApplyBQSR_connor:
    input:
        bam = "Connor/bam/{sample}.connor.dedup.bam",
        bai = "Connor/bam/{sample}.connor.dedup.bam.bai",
        recal = "Connor/quality_control/{sample}_pre_recalibration.grp",
    output:
        bam = "Connor/bam/{sample}.connor.recalibrated.bam",
        bai = "Connor/bam/{sample}.connor.recalibrated.bai",
    shell:
        """
        gatk --java-options -Xmx{mem}G ApplyBQSR \
        -R {ref} \
        -I {input.bam} \
        --bqsr {input.recal} \
        -O {output.bam}
        """
#-L {input.bed}
#bed = {bed}

rule post_recalibration_table_connor:
    input:
        bam = "Connor/bam/{sample}.connor.recalibrated.bam",
        bai = "Connor/bam/{sample}.connor.recalibrated.bai",
    output:
        "Connor/quality_control/{sample}_post_recalibration.grp"
    shell:
        """
        gatk --java-options -Xmx{mem}G BaseRecalibrator \
        -R {ref} \
        -I {input.bam} \
        --known-sites {dbsnp} \
        --known-sites {mills_1000g} \
        --known-sites {phase1_1000g} \
        -O {output}
        """

#-L {input.bed}
#bed = {bed}

rule CollectHsMetrics_connor:
    input:
        bam = "Connor/bam/{sample}.connor.recalibrated.bam",
        bai = "Connor/bam/{sample}.connor.recalibrated.bai",
        interval = {interval_list}
    output:
        "Connor/quality_control/{sample}_capture-targets_HsMetrics.txt"
    shell:
        """
        gatk --java-options -Xmx{mem}G CollectHsMetrics \
        --INPUT {input.bam} \
        --REFERENCE_SEQUENCE {ref} \
        --OUTPUT {output} \
        --BAIT_INTERVALS {input.interval} \
        --TARGET_INTERVALS {input.interval}
        """

rule CollectInsertSizeMetrics_connor:
    input:
        bam = "Connor/bam/{sample}.connor.recalibrated.bam"
    output:
        hist = "Connor/quality_control/{sample}.insert_size_histogram.pdf",
        metrics = "Connor/quality_control/{sample}.insert_size_metrics.txt"
    shell:
        """
        gatk --java-options -Xmx{mem}G CollectInsertSizeMetrics \
        -I {input.bam} \
        -H {output.hist} \
        -O {output.metrics}
        """


rule qualimap_connor:
    input:
        "Connor/bam/{sample}.connor.recalibrated.bam"
    output:
        html = "Connor/quality_control/{sample}/qualimapReport.html"
    params:
        outdir = directory("Connor/quality_control/{sample}")
    resources: cpus=4, mem=30000
    shell:
        """
        qualimap bamqc \
        -bam {input} \
        -nt 24 \
        -c \
        -sd \
        -gff {bed} \
        -outdir {params.outdir} \
        --java-mem-size=40G
        """

rule multiqc_connor:
    input:
        expand("quality_control/fastqc/{sample}_R1_fastqc.html", sample=SAMPLES),
        expand("quality_control/fastqc/{sample}_R1_fastqc.zip", sample=SAMPLES),
        expand("quality_control/fastqc/{sample}_R2_fastqc.html", sample=SAMPLES),
        expand("quality_control/fastqc/{sample}_R2_fastqc.zip", sample=SAMPLES),
        expand("Connor/quality_control/{sample}_pre_recalibration.grp", sample=SAMPLES),
        expand("Connor/quality_control/{sample}/qualimapReport.html", sample=SAMPLES),
        expand("Connor/quality_control/{sample}_capture-targets_HsMetrics.txt", sample=SAMPLES),
        expand("Connor/quality_control/{sample}_post_recalibration.grp", sample=SAMPLES)
    output:
        html = "Connor/quality_control/{famname}_connor_multiqc_report.html"
    params:
        indir1 = directory("Connor/quality_control"),
        indir2 = directory("quality_control/fastqc/"),
        config = resource_path + "multiqc_config.yaml"
    shell:
        """
        multiqc {params.indir1} {params.indir2} \
        -n {output.html} \
        -c {params.config}
        """


##################################
######       Version Log      ####
##################################

rule version_log:
    input:
        # markdup = expand("MarkDuplicates/quality_control/{famname}_multiqc_report.html", famname=FAMNAME),
        connor = expand("Connor/quality_control/{famname}_connor_multiqc_report.html", famname=FAMNAME)
    output:
        log="logs/exom_align_tagseq_combined_hg38-version-log_"+timeFormat+".txt"
    shell:
        """
        echo 'Script: exom_align_tagseq_combined_hg38' > {output};
        conda list --export >> {output}
        """
