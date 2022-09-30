#===================================================================================
#
# FILE: outlier_suppression.R
# USAGE: Apply outlier suppression to samples after they are split.
#       This first determines a ctDNA level across the loci, using which
#       outliers can be identified
#
#===================================================================================

## libraries
library(dplyr)
library(data.table)

# Read in command line variables
args = commandArgs(TRUE)

print("loading args")
data.path <- args[1]
file_prefix <- args[2]
code_dir <- args[3]

# Read in TAPAS functions file
print("source functions")
source(paste0("config/config2.R"))
source(paste0(code_dir, "/R/TAPAS_functions.R"))
source(paste0(code_dir, "/R/detection_functions.R"))


print(paste("loading", data.path))

load(data.path)
raw <-repolish(curr, outlier_suppression_threshold_settings)

print("saving")
save(raw, file = paste0("output_R/os/", file_prefix, ".os.Rdata"))

## only save files with mut reads
if (sum(raw$mutant) > 0){
  print("mutant reads present")
  nonzero <-filter(raw, mutant == 1)
  nonzero <- nonzero[!duplicated(nonzero$uniq),]

  save(nonzero, file = paste0("output_R/os_nonzero/", file_prefix, ".os.nonzero.Rdata"))
}
