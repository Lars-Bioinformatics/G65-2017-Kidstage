# WORK="/work/sduvarcall/G65-2017-Kidstage/Connor/cgp_bam2"
WORK="/work/sduvarcall/G65-2017-Kidstage/MarkDuplicates/cgp_bam"

SAMPLES, = glob_wildcards(WORK+"{sample}.bam")

# Resources - paths inside docker
# ref_fai = "/work/sduvarcall/bwa-0.7.13/reference/human_g1k_v37_decoy.fasta.fai"
ref_fai = "/work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.fasta.fai"

# onstart:
#     shell("mkdir -p " + OUTPUT)

rule all:
    input:
        expand(WORK+"{sample}.bam.bas", sample=SAMPLES),
        # expand(WORK+"{sample}.bam.bai", sample=SAMPLES)


rule bam_stats:
    input:
        bam=WORK+"{sample}.bam"
    output:
        bas=WORK+"{sample}.bam.bas"
    shell:
        """
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume=/work:/work \
        cgpwgs \
        bam_stats \
        -i {input} \
        -o {output} \
        -r {ref_fai} \
        -@ 12
        """

rule rename_bam_index:
    input:
        bai=WORK+"{sample}.bai"
    output:
        bai=WORK+"{sample}.bam.bai"
    shell:
        """
        cp {input} {output}
        """