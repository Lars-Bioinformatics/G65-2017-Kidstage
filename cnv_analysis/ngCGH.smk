__title__ = "Script for running ngCGH CNV analysis"
__author__ = "Lars Andersen <larsmew@gmail.com>"
__date__ = "19/01/2021"
__version__ = "1.0"

WORK = "/work/G65-2017-Kidstage/"

## Run from folder with recalibrated bam files ##
configfile: WORK+"samples-pairwise.yaml"

onstart:
    shell("mkdir -p "+WORK+"ngcgh")

# regions = "/work/sdukoldby/resources/hg38/MedExome_hg38_capture_targets.bed"
window_size = 1000

rule all:
    input:
        [expand(WORK+"ngcgh/{tumor}_vs_{normal}_tagseq-medexome_ngcgh_w{window_size}.txt",
                    window_size=window_size,
                    normal=config[pair]["normal"],
                    tumor=config[pair]["tumor"]) for pair in config]

rule ngcgh:
    input:
        normal = WORK+"Connor/bam/{normal}.connor.recalibrated.bam",
        tumor  = WORK+"Connor/bam/{tumor}.connor.recalibrated.bam",
    output:
        out = WORK+"ngcgh/{tumor}_vs_{normal}_tagseq-medexome_ngcgh_w{window_size}.txt"
    conda: "envs/ngcgh.yaml"
    threads: 1
    shell:
        """
        ngCGH {input.normal} {input.tumor} \
        -w {window_size} \
        -o {output}
        """
