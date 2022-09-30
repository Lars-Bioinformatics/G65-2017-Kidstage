
#===================================================================================
#
# FILE: INVAR4.R
# USAGE: run after INVAR3 to annotate on target data with SLX table information.
#        SLX table is a spreadsheet which contains sequencing run information.
#        An example is given in the example_files folder. Please see the wiki
#
#===================================================================================

library(plyr)
library(data.table)
library(dplyr)
library(ggplot2)
library(plotROC)

source("~/tapas_pipeline/R/TAPAS_functions.R")
source("config/config2.R")

## set variables
args = commandArgs(TRUE)

tumour_mutations.path <- args[1]
output_dir <- "output_R/"
input_dir <- "output_R/"
SLX_layout_path <- args[2] 

# load and prepare SLX layout table -----------------------------
SLX_layout <- read.csv(file = SLX_layout_path, header = TRUE)
SLX_layout$barcode. <- gsub(SLX_layout$barcode., pattern = "-", replacement = "_")
SLX_layout$SLX_barcode <- paste0(SLX_layout$SLX_ID, "_", SLX_layout$barcode.)

# read in on_target dataframe
study <- unlist(lapply(strsplit(list.files("output_R", pattern = "on_target.Rdata"), split = "\\."), "[[", 1))
setting <- gsub(list.files("output_R", pattern = "on_target.Rdata"), pattern = ".*(f0.._s[0-9].BQ_...MQ_..).*|raw", replacement = "\\1")

print(paste("loading", list.files("output_R", pattern = "on_target.Rdata")))
print(paste("setting = ", setting))

# load data_dt.on_target
load(list.files("output_R", pattern = "on_target.Rdata", full.names = T))

# read in tumour mutations
tumour_mutations <- read.csv(file = tumour_mutations.path)
error_file.path <- list.files("output_R", pattern = "FALSE.error_rates", full.names = T)
error_file.path <- error_file.path[!grepl(error_file.path, pattern = "locus")]

# check all patients are represented in tumour mutations file
print("patients NOT represented in tumour mutations file:")
print(unique(filter(SLX_layout, Study == study)$Patient)[!(unique(filter(SLX_layout, Study == study)$Patient) %in% unique(tumour_mutations$Patient))])

# check that all mutations are <100% allele fraction
if (all(tumour_mutations$tumour_AF <=1) == FALSE){
  print("error - 1 or more mutation loci have an allele fraction of >100%, please set to between 0 and 1.00 then try again")
  return(NULL)
}

print("loaded files successfully")
print(paste("proportion of controls = ", proportion_of_controls))
print(paste("max_background_mean_AF = ", max_background_mean_AF))

if( exists("outlier_suppression_threshold_settings") == FALSE){
  print("setting default OS setting as 0.05")
  outlier_suppression_threshold_settings <- 0.05
}

for (outlier_suppression_threshold in outlier_suppression_threshold_settings){
  print(paste("outlier suppression threshold = ", outlier_suppression_threshold))

  combined.polished <- TAPAS4(data_dt.on_target,
                        cosmic_threshold = cosmic_threshold,
                      proportion_of_controls = proportion_of_controls, # filter loci that are positive in >30% of control samples (need sufficient controls)
                      max_background_mean_AF = max_background_mean_AF,
                      outlier_suppression_threshold = outlier_suppression_threshold,
                      AF_threshold = 0.01, # do not use samples >5% AF for nonptspec background - contamination Risk
                      tumour_mutations, # tumour mutation dataframe
                      error_file = error_file.path, # error rates file from path
                      SLX_layout,
                      max_DP = 2000) # SLX layout table for sequencing annotation info

  # handle multiple patient-specific panel data; in this case it is only AVASTM
  if (grepl(study, pattern = "AVAST")== TRUE){
    print("handling multipanel data")
    combined.polished<- ensure_on_panel(combined.polished, tumour_mutations)
  }

  # save
  print("saving RDS")
  saveRDS(combined.polished, file = paste0("output_R/", study,".", setting, ".combined.rds"))
}
print("TAPAS4 complete")
