#===================================================================================
#
# FILE: size_annotation.R
# USAGE: annotate df with sizes, save as individual Rdata files for Eyal to run detection
#       For interactive mode: srun --mem=60000 --pty /home/wan01/software/R-3.3.2/bin/R --vanilla
#
#===================================================================================

# Read in command line variables
args = commandArgs(TRUE)

combined_polished.path <- args[1] #"/Users/wan01/shell_scripts/cluster1/TAPAS"  #
size_ann.file <- args[2] #'f0.9_s2.BQ_30.MQ_60'
final_prefix <- args[3]
code_dir <- args[4]
outlier_suppression <- args[5]

# echo variables
print(combined_polished.path)
print(size_ann.file)
print(final_prefix)
print(code_dir)
print(outlier_suppression)

# Read in TAPAS functions file
source(paste0(code_dir, "/R/TAPAS_functions.R"))
source(paste0("config/config2.R")) # note that this sets output_dir

#update output_dir
output_dir <- paste0("output_R/combined.polished.size")

# load combined_polished df
print(paste("loading", combined_polished.path))
combined.polished <- readRDS(combined_polished.path)

print(paste("working on ", size_ann.file))
combined.polished.size_ann <- equalise_size_counts(size_ann.csv = size_ann.file, combined.polished = combined.polished)

curr_SLX_barcode <- unique(combined.polished.size_ann$SLX_barcode)

# split by sample_name and ptspec/nonptspec
samples <- unique(combined.polished.size_ann$sample_name)

for (x in 1:length(samples)){
  print(paste("sample = ", samples[x]))
  curr <- filter(combined.polished.size_ann, sample_name == samples[x])
  ptspec_status <- unique(curr$data)

  print(paste("ptspec status: ", ptspec_status))

  # Run a check for both_strands to ensure that both strands really are represented
  curr <- annotate_with_both_strands(curr)
  print(paste("saving", samples[x], "data ==", ptspec_status))

  save(curr, file = paste0(output_dir, "/", final_prefix, ".", curr_SLX_barcode, ".", x, ".", ptspec_status, ".", outlier_suppression, ".combined.polished.size.Rdata"))
}
