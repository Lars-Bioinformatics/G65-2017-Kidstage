
#===================================================================================
#
# FILE: INVAR3.R
# USAGE: Set variables below, then run on error-suppressed INVAR BAM data that has
#       been processed using the make_mpileup.sh script
#       For interactive mode: srun --mem=60000 --pty /home/wan01/software/R-3.3.2/bin/R --vanilla
#
# DESCRIPTION: Raw INVAR data is parsed, filtered with 3 separate filters (MQSB, multiallelic, and strand bias).
#               Files are then split into target vs. nontarget based on loci of interest
#
#===================================================================================

# Read in command line variables
args = commandArgs(TRUE)

code_dir <- args[1] #"/Users/wan01/shell_scripts/cluster1/INVAR"  #
TAPAS.setting <- args[2] #'f0.9_s2.BQ_30.MQ_60'
mutations_bed  <- args[3] #MELR: "~/bed/M2/JXP0128_4_MELR_TAPAS.bed" #AVASTM: A1-7_panels.combined.tidy.180109.bed #LUCID:
SLX_layout_path <- args[4]
output_dir <- "output_R/"
off_target.split_output_dir <- "off_target.split" #dir for split off target output (split by SLX_barcode)

# Read in TAPAS functions file
source(paste0(code_dir, "/R/TAPAS_functions.R"))

# Read in filtering criteria
source(paste0("config/config2.R"))

load_existing_Rdata <- F

if(load_existing_Rdata == F){
  # Read in combined.txt.gz file
  print("reading in table from output_gz")
  raw <- read.table(paste0(gz_dir, TAPAS.setting, TAPAS.file_suffix), col.names = TAPAS.tsv_colnames)
  data <- raw

  # parse columns
  print("parsing columns to ensure the correct structure")
  data$REF <- as.character(data$REF)
  data$ALT <- as.character(data$ALT)
  data$MQSB <- as.numeric(as.character(data$MQSB))

  print(paste("summary of MQSB", summary(data$MQSB)))

  # Add columns
  data$AF <- (data$ALT_F + data$ALT_R) / data$DP
  data$SLX_barcode <- sprintf("%s_%s", data$SLX, data$BARCODE)
  data$uniq_pos <- sprintf("%s:%s", data$CHROM, data$POS)
  data$uniq <- sprintf("%s_%s", data$uniq_pos, data$SLX_barcode)

  # Annotate
  data$COSMIC <- data$COSMIC_MUTATIONS > cosmic_threshold
  data$SNP <- data[,grepl(colnames(data), pattern = "KG_AF")] > 0

  # annotate with on/off target
  print("reading in patient-specific loci from the BED")
  patient_specific <- read.table(mutations_bed, fill = T)
  patient_specific$uniq_pos <- paste0(patient_specific[,1], ":", patient_specific[,3])
  data$on_target <- data$uniq_pos %in% patient_specific$uniq_pos

  # remove soft-masked repeats, identified by lowercase
  print("Skipping loci that were soft-masked repeats")
  data <- data[!grepl(data$REF, pattern = "a|c|g|t"), ]
  data <- data[!grepl(data$ALT, pattern = "a|c|g|t"), ]

  # save on.target pre-filtering
  data_dt.prefilter <- data.table(data)

  # Get the on-target bases (containing patient-specific and non-patient-specific)
  print("Filter data data.table to collect on_target data, SNP == FALSE")
  data_dt.prefilter <- data_dt.prefilter[on_target == TRUE &
                     # prefilter data frame only filters out INDELS and off target reads.
                     #SNP == FALSE &
                     nchar(ALT) == 1 &
                     nchar(REF) == 1]

  print("saving data_dt.prefilter")
  save(data_dt.prefilter, file = paste0(output_dir, TAPAS.setting, ".on_target.prefilter.Rdata"))

  # Apply filters to blacklist loci
  data <- filter(data, DP < max_DP,
                       SNP == FALSE,
                       REF_R + REF_F >= min_ref_DP,
                       (COSMIC == FALSE & on_target == FALSE | on_target == TRUE))
  data <- blacklist.MQSB(data, individual_MQSB_threshold)
  data <- blacklist.multiallelic(data, n_alt_alleles_threshold, minor_alt_allele_threshold)
  print("Strand bias filter has been removed and will no longer be used.")


  # save data as an Rdata file
  dir.create(output_dir)
  print("saving data prior to filtering as Rdata file")
  save(data, file = paste0(output_dir, TAPAS.setting, ".Rdata"), compress = F)

} else{
  load(paste0(output_dir, TAPAS.setting, ".Rdata"))
}

## save on/off target bases as separate files
data_dt.raw <- data.table(data)

# Get the on-target bases (containing patient-specific and non-patient-specific)
data_dt.on_target <- data_dt.raw[on_target == TRUE &
                   SNP == FALSE &
                   nchar(ALT) == 1 &
                   nchar(REF) == 1]

save(data_dt.on_target, file = paste0(output_dir, TAPAS.setting, ".on_target.Rdata"))

#take only off_target bases, exclude INDELs
for (cosmic in c(FALSE, TRUE)){
  print(paste("COSMIC = ", cosmic))

  if (cosmic == FALSE){
    print("removing cosmic loci")
    data_dt.off_target <- data_dt.raw[on_target == FALSE &
                       SNP == FALSE &
                       COSMIC == FALSE &
                       nchar(ALT) == 1 &
                       nchar(REF) == 1]
     } else{
       print("NOT removing cosmic loci")
       data_dt.off_target <- data_dt.raw[on_target == FALSE &
                          SNP == FALSE &
                          nchar(ALT) == 1 &
                          nchar(REF) == 1]
     }

  # error rate per locus filters
  background_error.raw <- data_dt.off_target[,.(mut_sum=sum(ALT_F)+sum(ALT_R),
                                  total_DP=sum(DP),
                                  background_AF=(sum(ALT_F)+sum(ALT_R))/sum(DP)),
                                  by=list(SLX, BARCODE, REF, ALT, TRINUCLEOTIDE)]
  print.error_free_positions(data_dt.off_target)

  data_dt.off_target <- annotate_with_locus_error_rate(combined.filtered = data_dt.off_target, #data_dt.off_target,
                                             is.blood_spot = is.blood_spot,
                                             on_target = FALSE)

  # annotate with both strand information
  data_dt.off_target <- both_strands(data_dt.off_target)

  # calculate error rate with differnt settings
  # locus noise pass
  background_error.locus_noise <- data_dt.off_target[LOCUS_NOISE.PASS == TRUE,.(mut_sum=sum(ALT_F)+sum(ALT_R),
                                  total_DP=sum(DP),
                                  background_AF=(sum(ALT_F)+sum(ALT_R))/sum(DP)),
                                  by=list(SLX, BARCODE, REF, ALT, TRINUCLEOTIDE)]
  print.error_free_positions(filter(data_dt.off_target, LOCUS_NOISE.PASS == TRUE))

  # both strands
  data_dt.off_target.both_strands <- filter(data_dt.off_target, BOTH_STRANDS == TRUE) %>% data.table()
  background_error.both_strands <- data_dt.off_target.both_strands[,.(mut_sum=sum(ALT_F)+sum(ALT_R),
                                  total_DP=sum(DP),
                                  background_AF=(sum(ALT_F)+sum(ALT_R))/sum(DP)),
                                  by=list(SLX, BARCODE, REF, ALT, TRINUCLEOTIDE)]
  print.error_free_positions(data_dt.off_target.both_strands)

  # locus noise AND both strands
  data_dt.off_target.locus_noise.both_strands <- filter(data_dt.off_target, LOCUS_NOISE.PASS == TRUE, BOTH_STRANDS == TRUE) %>% data.table()
  background_error.locus_noise.both_strands <- data_dt.off_target.locus_noise.both_strands[,.(mut_sum=sum(ALT_F)+sum(ALT_R),
                                  total_DP=sum(DP),
                                  background_AF=(sum(ALT_F)+sum(ALT_R))/sum(DP)),
                                  by=list(SLX, BARCODE, REF, ALT, TRINUCLEOTIDE)]
  print.error_free_positions(data_dt.off_target.locus_noise.both_strands)

  error_rates <- list(pre_filter = background_error.raw,
      locus_noise_filter_only = background_error.locus_noise,
       both_strands_only = background_error.both_strands,
       locus_noise.both_strands = background_error.locus_noise.both_strands)

  save(error_rates, file = paste0("output_R/", TAPAS.setting, ".cosmic_", cosmic, ".error_rates.Rdata"))

  if(cosmic == TRUE){
    print("cosmic == true, saving off target files for ")
    ## save multiple raw Rdata objects of each off-target sample
    # for eyal
    dir.create(paste0(output_dir, off_target.split_output_dir))

    plyr::ddply(data_dt.off_target.locus_noise.both_strands, "SLX_barcode", function(x){
      curr_barcode <- unique(x$SLX_barcode)
      print(paste0("saving off-target, locus_noise.both_strands background data for: ", curr_barcode))

      save(x, file = paste0(output_dir,
                            off_target.split_output_dir, "/",
                            TAPAS.setting,
                            ".off_target.locus_noise.both_strands.",
                            curr_barcode,
                            ".Rdata"))
    })
  }
}

# save a dataframe of the filterted loci
print("saving a dataframe of the loci that were filtered out based on MQSB, locus noise etc.")
filtered_out <- data_dt.prefilter[!(data_dt.prefilter$uniq_pos %in% data$uniq_pos), ]
save(filtered_out, file = paste0("output_R/", TAPAS.setting, ".on_target.filtered_out.Rdata"))
