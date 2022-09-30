#===================================================================================
#
# FILE: size_annotation3.R
# USAGE: generate size characterisation data frame. This only takes the loci that have passed all filters
#         The resulting dataframe is used in the GLRT code
#
#===================================================================================

# Read in command line variables
args = commandArgs(TRUE)

final_prefix <- args[1]
SLX_table.path <- args[2]
code_dir <- args[3]

# Read in TAPAS functions file
source(paste0(code_dir, "/R/TAPAS_functions.R"))
source(paste0("config/config2.R")) # note that this sets output_dir

#update output_dir
output_dir <- paste0("output_R")

# read SLX table
SLX_layout <- read.csv(file = SLX_table.path,
    header = TRUE)
SLX_layout$barcode. <- gsub(SLX_layout$barcode., pattern = "-", replacement = "_")
SLX_layout$SLX_barcode <- paste0(SLX_layout$SLX_ID, "_", SLX_layout$barcode.)

files <- list.files(path = "output_R/os", pattern = "Rdata", full.names = T)

# load in files
output <- data.frame()

for (x in 1:length(files)){
  print(paste("running file:", x, "out of", length(files)))
  load(files[x])
  curr_dt <- data.table(raw)

  curr_dt.filtered <- curr_dt[PASS == TRUE & BOTH_STRANDS == TRUE & LOCUS_NOISE.PASS == TRUE ][,c("SLX_barcode", "sample_type", "sample_name", "data", "mutant", "size")]

  curr_dt.filtered$mutant <- curr_dt.filtered$mutant == 1

  print("summarising")
  summary <- curr_dt.filtered[,.(total = .N),
                                by = list(SLX_barcode, sample_name, sample_type, data, size, mutant)]

  output <- rbind(summary, output)
}

output <- left_join(output, SLX_layout[,c('case_or_control', 'SLX_barcode')], by = "SLX_barcode") %>% data.table()

print("saving raw size characterisation file")
saveRDS(output, paste0(output_dir, "/", final_prefix, ".raw.size_characterisation.rds"))

## summarise across samples
summary <- output[,.(count = sum(total)),
                                by = list(data, case_or_control, mutant, size)]

totals <- plyr::ddply(summary, c("data", "case_or_control", "mutant"), function(x){
                      data.frame(total = sum(x$count)) })


summary <- left_join(summary, totals, by = c("data", "case_or_control", "mutant"))

summary$proportion <- summary$count / summary$total
colnames(summary)[3] <- "mut"

print("saving summarised size characterisation file")
saveRDS(summary, paste0(output_dir, "/", final_prefix, ".size_characterisation.rds"))
