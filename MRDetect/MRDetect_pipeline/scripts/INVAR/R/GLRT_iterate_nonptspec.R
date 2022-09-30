#===================================================================================
#
# FILE: GLRT_iterate_nonptspec.R
# USAGE: GLRT code, run in parallel by using GLRT.sh. This code iterates the nonptspec
#        samples when calculating GLRT in order to more accurately determine the distribution of
#         INVAR scores in the controls
#
#===================================================================================

## libraries
library(dplyr)
library(data.table)

# Read in command line variables
args = commandArgs(TRUE)

print("loading args")
data.path <- args[1]
size_data.path <- args[2]
detection_functions <- args[3]
smooth <- args[4]
size_data.path.prefix <- args[5]
final_prefix <- args[6]

# determine the outlier suppression level
outlier_suppression <- gsub(data.path, pattern =".*\\.(0\\.0[0-9]+)\\..*", replacement = "\\1")
print(paste("outlier suppression = ", outlier_suppression))

# Read in TAPAS functions file
print("sourcing detection functions")
source(detection_functions)
source(paste0("config/config2.R"))
source("~/tapas_pipeline/R/TAPAS_functions.R")

# check for variables
if(exists("downsample_loci") ==FALSE){
  print("not downsampling loci. [default]")
  downsample_loci <- NA
} else if (exists("downsample_loci") ==TRUE){
  print("downsample loci exists. downsampling")
} else{
  print("error")
  return(NULL)
}

print("loading in files")
load(as.character(data.path))
size_characterisation <- readRDS(size_data.path) %>% data.table()

## select columns depending on what data has been used
if ("data" %in% colnames(size_characterisation) == TRUE){
  print("size characterisation dataframe has additonal columns, filtering")
  size_characterisation <- filter(size_characterisation, data == "ptspec") %>%
                              select(mut, total, size, count, proportion)
}


# set min and max fragment length
# N.B. this acts as a size selection
min_length <- 60
max_length <- 300

if (exists("curr") == TRUE){
  print("running regular GLRT")
  raw <- curr
} else{
  print("running repolished")

}

if (all(grepl(raw$sample_name, pattern = "MIXTURE")) == FALSE |
          (unique(raw$data) == "nonptspec")){
  print("not downsampling loci, default")
  downsample_loci <- NA
}

### calculate IMAFv2 (BS_uAF; trinucleotide background-subtracted mean AF)
IMAFv2 <- calculate_IMAFv2(raw, is.blood_spot)


for (PASS_boolean in c(TRUE, FALSE)){
    for (BOTH_STRANDS_boolean in c(TRUE)){
        for (LOCUS_NOISE.PASS_boolean in c(TRUE)){

        # get raw filtered object again
      raw.filtered <- raw
      print(paste("current condition (PASS; BOTH_STRANDS; LOCUS_NOISE.PASS):", PASS_boolean, BOTH_STRANDS_boolean, LOCUS_NOISE.PASS_boolean))

      print(paste("number of rows of curr:", nrow(raw.filtered)))

      #filters in a stepwise fashion, only filtering if the condition == TRUE
      if(PASS_boolean == TRUE ){
        print("filtering for PASS")
        raw.filtered <- filter(raw.filtered, PASS == PASS_boolean)
      } else{
        print("not filtering for repolished loci here as they will be repolished later")
      }

      #print(paste("number of rows of curr:", nrow(curr)))
      if(BOTH_STRANDS_boolean == TRUE){
        print("filtering for BOTH STRANDS")
        raw.filtered <- filter(raw.filtered, BOTH_STRANDS == BOTH_STRANDS_boolean)
      }

      #print(paste("number of rows of curr:", nrow(curr)))
      if(LOCUS_NOISE.PASS_boolean == TRUE & is.na(downsample_loci[1])){
        print("filtering for LOCUS NOISE")
        raw.filtered <- filter(raw.filtered, LOCUS_NOISE.PASS == LOCUS_NOISE.PASS_boolean)
      }

      # leave-one-out style exclusion of mutant reads from this particular sample from the
      # large size_characterisation df to ensure that approach is non-circular
      print("starting leave_one_out code")
      size_characterisation <- leave_one_out.size(raw.filtered = raw.filtered,
                   size_characterisation = size_characterisation)


      # tumour AF and size filters
      data <- filter(raw.filtered,
                          tumour_AF > 0,
                         size > min_length,
                         size <= max_length)

     #print(paste("number of rows of data:", nrow(curr)))
      # for error classes where there are zero mutant reads, add one just so that the background error estimate isn't zero
      if (length(which((data$background_AF == 0))) > 0){
        print("some error classes have zero error rate, increase error rate to 1/background.DP error rate")
        print(unique(data[data$background_AF == 0, 'mut_class']))
        data[data$background_AF == 0, 'background_AF'] <- 1/data[data$background_AF == 0, 'background.DP']
      }


    data.pre_downsample <- data

    # This loop downsamples loci for the mixture experiment. By default, downsample loci is set to NA
    # only set downsample loci variable if you want to downsample loci
    for (x in downsample_loci){

      # In silico downsampling experiment of LOCI
      if( all(!is.na(downsample_loci)) & all(grepl(data.pre_downsample$sample_name, pattern = "MIXTURE"))) {
        print("JW downsample experiment - downsampling loci")
        data <- downsample_loci_fn(data.pre_downsample)
      }

      # determine whether there are any mutant reads in the dataframe for calculation of ctDNA
      if (nrow(filter(data, mutant == 1)) > 0){
        mutant_reads_present <- T
      } else{
        mutant_reads_present <- F
      }

      # only count mutant reads once for accurate ctDNA quantification
      data.ordered <- data[order(data$uniq, data$size),]
      data1 <- data.ordered[seq(1, nrow(data.ordered), by = 2),]

      if (unique(data1$data) == "nonptspec"){
        n_iterations <- 10
      } else{
        n_iterations <- 1
      }
        for (y in 1:n_iterations){

          if (unique(data1$data) == "nonptspec" & y > 1){
            print("resampling (iteration > 1)")
            data1 <- data1[sample(nrow(data1), nrow(data1), replace = T),]
          } else {
            print ("not resampling, iteration == 1")
          }

            print(paste("Run characterise likelihood ratio for sample function"))
            output <- calculate_likelihood_ratio_for_sample(data1, size_characterisation, min_length, max_length, use_size = TRUE, smooth = smooth, size_data.path.prefix, final_prefix, only_weigh_mutants)
            output.no_size <- calculate_likelihood_ratio_for_sample(data1, size_characterisation, min_length, max_length, use_size = FALSE, smooth = smooth, size_data.path.prefix, final_prefix, only_weigh_mutants)

            # calculate ctDNA level
            print(paste("INVAR_SCORE.using_size,DP,AF_p,IMAF,Patient,pt_mutation_belongs_to,PASS,BOTH_STRANDS,LOCUS_NOISE.PASS,smooth,size_characterisation.path,CONTAMINATION_RISK.PASS,outlier_suppression,null_likelihood,alt_likelihood,mutant_reads_present,downsampled,mut_sum,iteration:", output[[1]], nrow(data), output[[2]], IMAFv2, unique(data$Patient), unique(data$pt_mutation_belongs_to), PASS_boolean, BOTH_STRANDS_boolean, LOCUS_NOISE.PASS_boolean, smooth, size_data.path, unique(data$CONTAMINATION_RISK.PASS), outlier_suppression, output[[3]], output[[4]], mutant_reads_present, x, sum(data$mutant),y))
            print(paste("INVAR_SCORE.no_size,DP,AF_p,IMAF,Patient,pt_mutation_belongs_to,PASS,BOTH_STRANDS,LOCUS_NOISE.PASS,smooth,size_characterisation.path,CONTAMINATION_RISK.PASS,outlier_suppression,null_likelihood,alt_likelihood,mutant_reads_present,downsampled,mut_sum,iteration:", output.no_size[[1]], nrow(data), output.no_size[[2]], IMAFv2, unique(data$Patient), unique(data$pt_mutation_belongs_to), PASS_boolean, BOTH_STRANDS_boolean, LOCUS_NOISE.PASS_boolean, smooth, size_data.path, unique(data$CONTAMINATION_RISK.PASS), outlier_suppression, output[[3]], output[[4]], mutant_reads_present, x, sum(data$mutant),y))

            rm(output)
            rm(output.no_size)
          }
      }
    }
  }
}
