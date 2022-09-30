__title__ = "Pipeline for Somatic Variant Calling and Copy Number Analysis"
__author__ = "Jakob Jersild Nielsen <jjn@inguz.dev>"
__date__ = "2020/11/09"
__version__ = "1.1"

import time, os, sys, glob, pandas as pd, multiprocessing as mp
from ruamel_yaml import YAML


#########################################################
####                      Configs                    ####
#########################################################
configfile: "config.yaml"

input_dir = config["inputdir"]
plasma_dir = config["plasmadir"]
output_dir = config["outputdir"]
output_raw = output_dir+"raw/"
output_variants_notfiltered = output_dir+"variants_notfiltered/"
output_variants_filtered = output_dir+"variants_filtered/"
output_bam_dir = input_dir+"repsmerged/"

ref = config["ref"]
gnomad = config["gnomad"]
medexome = config["medexome"]
dbsnp = config["dbsnp"]
common_variants = config["common_variants"]
centromeres = config["centromeres"]
pon_1000g = config["pon_1000g"]
mem = config["mem"]
threadsparallel = config["threadsparallel"]
bamsuffix = config["samplesubset"]["bamsuffix"]
target_regions = config["target_regions"]
INTERVALS, = glob_wildcards(target_regions+"{interval}.interval_list")
INTERVALS = sorted(INTERVALS)

TYPEVAR = ["snvs","indels"]

VARIANTCALLEROUTPUTS = []
VARIANTCALLERS = []
yaml = YAML().load(open("config.yaml"))["variantcaller"]
for key, value in yaml.items():
    if config["variantcaller"][str(key)]["include"]:
        VARIANTCALLERS.extend([str(key)])
        if config["variantcaller"][str(key)]["variants"] == "snvs":
            VARIANTCALLEROUTPUTS.extend([str(key)+"_snvs.vcf"])
        elif config["variantcaller"][str(key)]["variants"] == "indels":
            VARIANTCALLEROUTPUTS.extend([str(key)+"_indels.vcf"])
        elif config["variantcaller"][str(key)]["variants"] == "both":
            VARIANTCALLEROUTPUTS.extend([str(key)+"_snvs.vcf", str(key)+"_indels.vcf"])
        elif config["variantcaller"][str(key)]["variants"] == "bothcombined":
            VARIANTCALLEROUTPUTS.extend([str(key)+"_bothcombined.vcf"])

SOMATICSEQ = ["Consensus.sSNV.vcf","Consensus.sINDEL.vcf","Ensemble.sSNV.tsv","Ensemble.sINDEL.tsv"]
SOMATICSEQ_MINCALLERS=[int(round(len(VARIANTCALLERS)/2)-1),int(round(len(VARIANTCALLERS)/2)),int(round(len(VARIANTCALLERS)/2)+1)]
# ENSEMBLE = ["somaticseq","somaticseq_mincallers"+str(SOMATICSEQ_MINCALLERS[0]),"somaticseq_mincallers"+str(SOMATICSEQ_MINCALLERS[1]),"somaticseq_mincallers"+str(SOMATICSEQ_MINCALLERS[2]),"mrdetect","mutect2","varscan2"]
# ENSEMBLE = ["somaticseq","mrdetect","mutect2","varscan2"]
ENSEMBLE = ["mrdetect"]
print(SOMATICSEQ_MINCALLERS)
print(ENSEMBLE)

if config["tsvOnStart"]:
    shell("conda run --name r Rscript code/samples.R")

MATCHED_TUMOR_PLASMA_NORMAL = pd.read_table(config["tumor-plasma-normal"], sep="\t", header=None)
MATCHED_TUMOR_NORMAL = pd.read_table(config["tumor-normal"], sep="\t", header=None)
# print(MATCHED_TUMOR_NORMAL)
ENSEMBLE_DATAFRAME = pd.read_table(config["ensemble"], sep="\t", header=None)
# print(ENSEMBLE_DATAFRAME)

def removeprefix(self: str, prefix: str, /) -> str:
    if self.startswith(prefix):
        return self[len(prefix):]
    else:
        return self[:]
def removesuffix(self: str, suffix: str, /) -> str:
    if self.endswith(suffix):
        return self[:-len(suffix)]
    else:
        return self[:]

# PLASMA_SAMPLES = ["G35-P"+str(i)+"-PLASMA_thruplex-tag-seq-hv_HVFFCDSXY" for i in range(1,11)]
# PLASMA_SAMPLES = []
# for SAMPLE in glob.glob(plasma_dir+"*.bam"):
#     PLASMA_SAMPLES.append(removesuffix(removeprefix(SAMPLE, plasma_dir), bamsuffix))
# print(PLASMA_SAMPLES)

ALL_PLASMA = [MATCHED_TUMOR_PLASMA_NORMAL[2][i] for i in MATCHED_TUMOR_PLASMA_NORMAL.index]
print(ALL_PLASMA)
print(len(ALL_PLASMA))

log_file = "log_file_somatic.txt"


#########################################################
####                       Setup                     ####
#########################################################
# Timing
totim = time.time()
timeFormat = "%Y_%m_%d:%X" # year, month, day, time H:M:S


#########################################################
####  Define workflow start, stop and error actions  ####
#########################################################
onstart:
    shell("mkdir -p "+output_dir)
    shell("mkdir -p "+output_raw)
    shell("mkdir -p "+output_variants_notfiltered)
    shell("mkdir -p "+output_variants_filtered)
    shell("mkdir -p "+output_variants_merged)
    shell("mkdir -p "+output_bam_dir)

#########################################################
####                  Run All Rules                  ####
#########################################################
rule all:
    input:
    #    variantcalleroutput=[expand(expand(output_variants_filtered+"{{tumor}}_vs_{{normal}}/{variantcalleroutput}", 
    #        variantcalleroutput=VARIANTCALLEROUTPUTS),
    #        normal=MATCHED_BIOPSY_NORMAL[1][i],
    #        tumor=MATCHED_BIOPSY_NORMAL[0][i]) for i in MATCHED_BIOPSY_NORMAL.index],
    #    ensemble=[expand(output_dir+"variants_final/{ensemble}/{tumor}_vs_{normal}/{somaticseq}",
    #        normal=MATCHED_BIOPSY_NORMAL[1][i],
    #        tumor=MATCHED_BIOPSY_NORMAL[0][i],
    #        somaticseq=SOMATICSEQ,
    #        ensemble=ENSEMBLE) for i in MATCHED_BIOPSY_NORMAL.index],
    #    somaticseq=[expand(expand(output_dir+"variants_final/somaticseq_mincallers{SomaticSeq_mincallers}/{{tumor}}_vs_{{normal}}/{{somaticseq}}",
    #        SomaticSeq_mincallers=SOMATICSEQ_MINCALLERS),
    #        normal=MATCHED_BIOPSY_NORMAL[1][i],
    #        tumor=MATCHED_BIOPSY_NORMAL[0][i],
    #        somaticseq=SOMATICSEQ) for i in MATCHED_BIOPSY_NORMAL.index],
        # varscancna=[expand(output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}_vs_{normal}.cnv.bed", 
        #     normal=MATCHED_TUMOR_NORMAL[1][i], 
        #     tumor=MATCHED_TUMOR_NORMAL[0][i]) for i in MATCHED_TUMOR_NORMAL.index],
        # vafvafcompare=[expand(output_dir+"vaf-vaf-plots/ensemble/{tumor}_vs_{normal}/complete",
        #     normal=MATCHED_BIOPSY_NORMAL[1][i],
        #     tumor=MATCHED_BIOPSY_NORMAL[0][i]) for i in MATCHED_BIOPSY_NORMAL.index],
        # mrdetectsnv=[expand(output_dir+"variants_final/{ensemble}/{tumor}_vs_{plasma}_vs_{normal}/MRDetectSNV_complete",
        mrdetectsnv=[expand(output_dir+"variants_final/{ensemble}/{tumor}_vs_{normal}/MRDetectSNV_complete",
            normal=MATCHED_TUMOR_PLASMA_NORMAL[1][i],
            tumor=MATCHED_TUMOR_PLASMA_NORMAL[0][i],
            # plasma=MATCHED_TUMOR_PLASMA_NORMAL[2][i],
            ensemble=ENSEMBLE) for i in MATCHED_TUMOR_PLASMA_NORMAL.index],
        # mrdetectcna_samplerunner=[expand(output_dir+"raw/varscan2cna/mrdetectcna/max150/{tumor}_vs_{normal}/{plasma}/{plasma}{suffix}",
        #     normal=MATCHED_TUMOR_NORMAL[1][i],
        #     tumor=MATCHED_TUMOR_NORMAL[0][i],
        #     plasma=PLASMA_SAMPLES,
        #     suffix=[".win500.delsdups.perwindow.median",".win500.neu.perwindow.median"]) for i in MATCHED_TUMOR_NORMAL.index],
#        meanvaf_bed=[expand(output_dir+"variants_final/all.bed")],
#        meanvaf_rc=output_dir+"meanvaf/MeanVAF_readcount_complete",
#        meanvaf_ann=[expand(output_dir+"meanvaf/{ensemble}/{tumor}_vs_{normal}/MeanVAF_annotator_complete",
#            normal=MATCHED_BIOPSY_NORMAL[1][i],
#            tumor=MATCHED_BIOPSY_NORMAL[0][i],
#            ensemble=ENSEMBLE) for i in MATCHED_BIOPSY_NORMAL.index],
#        meanvafsnvplots=[expand(output_dir+"plots/meanVAF-Z-score-barplots/snv_{ensemble}.pdf",
#            ensemble=ENSEMBLE)],
#        binaryvafsnvplots=[expand(output_dir+"plots/BinaryVAF-Z-score-barplots/snv_{ensemble}.pdf",
#            ensemble=ENSEMBLE)],
#        mrdetectsnvplots=[expand(output_dir+"plots/MRDetect-Z-score-barplots/snv_{ensemble}.pdf",
#            ensemble=ENSEMBLE)],
    # resources:
    #     mem_mb=12000
    # threads: 
    #     1

#########################################################
####                   Load Rules                    ####
#########################################################

# Variant Callers
for variantcaller in VARIANTCALLERS:
    if config["variantcaller"][variantcaller]["include"]:
        include: "rules/variantcallers/"+variantcaller+".smk"

# Filter Calls
include: "rules/filtercalls.smk"

# Create ensembles of variants
include: "rules/somaticseq.smk"

# Copy Number Analyses
include: "rules/cna/varscan2.smk"

# VAF-VAF plots
#include: "rules/vafvafplots.smk"

# MRDetect
include: "rules/mrdetect/mrdetect.smk"
include: "rules/mrdetect/mrdetectcna.smk"
include: "rules/mrdetect/snv-barplots.smk"

# meanVAF
include: "rules/meanVAF/meanVAF.smk"
include: "rules/meanVAF/snv-barplots.smk"
