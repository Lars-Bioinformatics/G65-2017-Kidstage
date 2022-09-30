#===================================================================================
#
# FILE: annotate_with_outlier_suppression.R
#
# DESCRIPTION: This script takes all of the nonzero lines from each of the split
#             data frames, so that mutant loci can be examined in one dataframe
#
#===================================================================================

## libraries
library(dplyr)
library(data.table)

# Read in command line variables
args = commandArgs(TRUE)
source("config/config2.R")

# get study and setting
study <- unlist(lapply(strsplit(list.files(path = "output_R", pattern = "on_target.prefilter.Rdata"), split = "\\."), "[[", 1))
setting <- gsub(list.files(path = "output_R", pattern = "on_target.prefilter.Rdata"), pattern = ".*(f0.._s[0-9].BQ_...MQ_..).*|raw", replacement = "\\1")

setwd("output_R/os_nonzero")

files <- list.files(pattern = "Rdata")

output <- data.frame()

for (x in 1:length(files)){
  load(files[x])
  print(paste(x, "of", length(files)))
  print(files[x])

  nonzero_fails <- select(nonzero, PASS, uniq) %>% distinct() %>% filter(PASS == FALSE)

  output <- rbind(nonzero_fails, output)
}

print(paste("the number of loci failing outlier suppression:", nrow(output)))
#saveRDS(output, file = paste0("../", study,".", setting, "combined.oulier_suppression", outlier_suppression_threshold, ".rds"))

# update combined RDS with PASS column based on 'output' dataframe
# check the number of files present
print("Reset working directory")
setwd("../../")

combined_rds.path <- list.files(path = "output_R", pattern= "combined", full.names=T)
print(combined_rds.path)
combined_rds.path <- combined_rds.path[grep(combined_rds.path, pattern = "rds")]

print("taking only the rds")
print(combined_rds.path)

if (length(combined_rds.path) ==1){
  print("path is correct, reading file...")
  combined.polished  <- readRDS(file = combined_rds.path)

  # ensure that we start with all loci passing
  combined.polished$PASS <- TRUE

  # annotate with loci that failed
  combined.polished[combined.polished$uniq %in% output$uniq, 'PASS'] <- F
  print("re-annotated combined.polished data.frame")

  print("tabling the PASS column:")
  print(table(combined.polished$PASS))
} else{
  print(combined_rds.path)
  print("There are multiple files (or none) matching combined*rds, please remove the unnecessary ones then retry")
  return(NULL)
}

## remove additional columns present in multi-panel data
combined.polished <-select(combined.polished, -starts_with("panel_number"))


# save the combined polished data.frame
# This data frame can be used to look at individual loci and their potential signal (annotated with outlier suppression)
print("saving RDS (now updated with a PASS column based on OS)")
print(paste0("Will use study = ", study, "Setting = ", setting, "and outlier threshold = ", outlier_suppression_threshold_settings))
saveRDS(combined.polished, file = paste0("output_R/", study,".", setting, ".combined.os_", outlier_suppression_threshold_settings ,".rds"))

print("done!")
