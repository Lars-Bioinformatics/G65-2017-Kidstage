#---
#title: "TAPAS pipeline functions"
#author: "JCM Wan and K Heider"
#---

## Load libraries --------------------------------------
if(!require(plyr) | !require(ggplot2) | !require(dplyr) | !require(readr) | !require(plotROC) | !require(data.table) | !require(IRanges) | !require(ggpubr)) {
  print("one or more packages are missing, installing them now..")
  install.packages(c("plyr", "dplyr", "ggplot2", "readr", "plotROC", "data.table", "ggpubr"))
}

## Filters for use on cluster  --------------------------------------
# MQSB
blacklist.MQSB <- function(data, individual_MQSB_threshold){
  #print(paste("starting rows (PRE-MQSB):", nrow(data)))
  print(paste("number of unique loci (pre MQSB)", length(unique(data$uniq_pos))) )

  print(paste("applying an invidiual MQSB threshold of:", individual_MQSB_threshold))
  data <- filter(data, MQSB > individual_MQSB_threshold)

  #print(paste("rows after filtering (POST-MQSB):", nrow(data)))
  print(paste("number of unique loci (post MQSB)", length(unique(data$uniq_pos))) )
  return(data)
}

# tidy multiallelic sites
blacklist.multiallelic <- function(data, n_alt_alleles_threshold = 3, minor_alt_allele_threshold = 2){
  print(paste("starting rows:", nrow(data)))
  nonzero <- filter(data, AF>0)

  # Determine the number of alt alleles per uniq_pos
  multiallelic <- plyr::ddply(nonzero, c("uniq_pos", "ALT"), function(x){
    mut_sum <- sum(x$ALT_F)+sum(x$ALT_R)
    data.frame(mut_sum)
  })

  # Identify loci with >1 alt alleles (e.g. both A>C and A>T at a position)
  multiallelic <- multiallelic[multiallelic$uniq_pos %in% multiallelic$uniq_pos[which(duplicated(multiallelic$uniq_pos))],]

  multiallelic2 <- plyr::ddply(multiallelic, "uniq_pos", function(x){
    data.frame(n_alt_alleles = nrow(x), min = min(x$mut_sum), max = max(x$mut_sum))
  })

  multiallelic_blacklist <- filter(multiallelic2, n_alt_alleles == n_alt_alleles_threshold,
                                   min >= minor_alt_allele_threshold,
                                   max >= minor_alt_allele_threshold)$uniq_pos

  print(paste("length of multiallelic blacklist =", length(multiallelic_blacklist)))

  # save the list of multiallelic loci for reference
  save(multiallelic_blacklist, file = "temp/multiallelic_blacklist.Rdata")

  data <- data[!(data$uniq_pos %in% multiallelic_blacklist),]

  print(paste("rows after filtering:", nrow(data)))
  return(data)
}

# strand bias
blacklist.strand_bias <- function(data, FR_ratio.min_reads = 10, FR_ratio.threshold = 0.2){
  print(paste("starting rows:", nrow(data)))

  print('calculating the sum of mut reads on F and R strand to assess F/R bias')
  nonzero <- filter(data, AF > 0)

  ## Assess the extent of strand bias
  ## Look for loci with a high ratio between the number of mutant reads on the F vs. R strand.
  sum_FR<- plyr::ddply(nonzero, 'uniq_pos', function(x){
    ALT_F <- sum(x$ALT_F)
    ALT_R <- sum(x$ALT_R)
    data.frame(ALT_F, ALT_R, bias = (ALT_F/ALT_R))
  })

    print('applying strand bias threshold')
    strand_bias  <- sum_FR[(sum_FR$bias > (1/FR_ratio.threshold) | sum_FR$bias < FR_ratio.threshold ) & (sum_FR$ALT_F + sum_FR$ALT_R >= FR_ratio.min_reads), 1]

  # save the list of strand bias loci for reference
  save(strand_bias, file = "temp/strand_bias_blacklist.Rdata")

  data <- data[!(data$uniq_pos %in% strand_bias), ]

  print(paste("rows after filtering:", nrow(data)))
  return(data)
}

## binom function for outlier suppression
binom <- function(x, p, n){
  if(x > 0){
    #print("non-zero p-value, running binom.test")
    return(binom.test(x, n, p, alternative="greater")$p.value)
  } else{
    return(1)
  }
}

repolish <- function(curr, os_threshold){
  print("Applying outlier suppression")
  print(paste("AF before: ", unique(curr$average)))
  raw <- curr

  # do not include loci with AF>0.25 or with >10mutant reads as we are trying to detect MRD. If there are loci with high AF, then there will also be some loci with lower AF (given sufficient loci).
  curr <- filter(curr, LOCUS_NOISE.PASS == TRUE, BOTH_STRANDS == TRUE,AF< 0.25, mut_sum < 10, tumour_AF > 0) #

  if (length(which((curr$background_AF == 0))) > 0){
    print("some error classes have zero error rate, increase error rate to 1/background.DP error rate")
    print(unique(curr[curr$background_AF == 0, 'mut_class']))
    curr[curr$background_AF == 0, 'background_AF'] <- 1/curr[curr$background_AF == 0, 'background.DP']
  }

  raw$average <- estimate_p_EM(M = curr$mutant,
                                R = rep(1, length(curr$mutant)),
                                AF = curr$tumour_AF,
                                e = curr$background_AF,
                                initial_p = 0.01,
                                iterations = 200)

 # estimate p for outlier suppression
 p_estimate <- max(unique(raw$average), weighted.mean(curr$AF, curr$tumour_AF))

 # if no mutant reads, AF estimate has to be zero
  if(sum(curr$mutant) == 0){
    p_estimate <- 0
  }

# now take the binom probability of seeing N mutant reads given p (ctDNA estimate in sample)
 raw$average <- p_estimate
 raw$prob <- mapply(binom, x = raw$ALT_F + raw$ALT_R,
                      n = raw$DP,
                      p = p_estimate)

  # calculate p-value threshold with a bonferroni correction for the number of loci tested
  p_threshold <- os_threshold/length(unique(curr$uniq_pos))
  raw$p_threshold <- p_threshold

  # samples pass the filter if they have p-values than that defined here.
  raw$PASS <- (raw$prob > p_threshold)

  return(raw)
}

# Annotate dataframe with whether the mutant read was observed in both F + R
both_strands <- function(INVAR_input){
  INVAR_input$BOTH_STRANDS <- (INVAR_input$ALT_R > 0 & INVAR_input$ALT_F > 0) | INVAR_input$AF == 0
  return(INVAR_input)
}

# count how many error free base positions there are in the data
print.error_free_positions <- function(data_dt.off_target.FILTERED){
  data_dt.off_target.FILTERED <- data.table(data_dt.off_target.FILTERED)
  error_free_positions <- (1-(nrow(data_dt.off_target.FILTERED[AF> 0])/nrow(data_dt.off_target.FILTERED)))*100
  print(paste("error-free positions:", error_free_positions))
}

# Annotate with locus error rate
# Locus error rate = overall background error rate per locus, aggregated across control samples
annotate_with_locus_error_rate <- function(combined.filtered,
                                           proportion_of_controls = 0.1,
                                           max_background_mean_AF = 0.01,
                                           is.blood_spot = F,
                                           on_target = TRUE, plot = T){

  if (is.blood_spot == TRUE){
    print("sWGS/blood spot data, do not set a max_background_mean_AF value as it is not appropriate in the low unique depth setting")
    max_background_mean_AF <- 1
  }

  if (on_target == TRUE){
    print("running code on on-target bases, only using nonptspec dataframe")
    patient_SLX_barcode_INVAR4 <- read.csv(file = SLX_layout_path, header = TRUE)
    patient_SLX_barcode_INVAR4 <- filter(patient_SLX_barcode_INVAR4, case_or_control == "case")
    patient_SLX_barcode_INVAR4$barcode. <- gsub(patient_SLX_barcode_INVAR4$barcode., pattern = "-", replacement = "_")
    patient_SLX_barcode_INVAR4$SLX_barcode <- paste0(patient_SLX_barcode_INVAR4$SLX_ID, "_", patient_SLX_barcode_INVAR4$barcode.)
    patient_SLX_barcode_INVAR4 <- patient_SLX_barcode_INVAR4[,"SLX_barcode"]
    error_rate.raw <<- filter(combined.filtered, data == "nonptspec") %>%
    filter(SLX_barcode %in% patient_SLX_barcode_INVAR4)%>%
      plyr::ddply(c("uniq_pos", "mut_class", "TRINUCLEOTIDE", "pt_mutation_belongs_to", "COSMIC"), function(x){
        mut_sum <- sum(x$ALT_F) +sum(x$ALT_R)
        DP <- sum(x$DP)
        n_samples <- length(unique(x$SLX_barcode))

        n_samples_with_signal <- filter(x, ALT_F + ALT_R > 0)$SLX_barcode %>% unique() %>% length()

        data.frame(background_AF = mut_sum / DP, mut_sum, DP, n_samples, n_samples_with_signal)


      })

      save(error_rate.raw, file = "output_R/locus_error_rates.on_target.Rdata")
  } else{
    print("running code on cluster, using all off-target bases")
    patient_SLX_barcode <- read.csv(file = SLX_layout_path, header = TRUE)
    patient_SLX_barcode <- filter(patient_SLX_barcode, case_or_control == "case")
    patient_SLX_barcode$barcode. <- gsub(patient_SLX_barcode$barcode., pattern = "-", replacement = "_")
    patient_SLX_barcode$SLX_barcode <- paste0(patient_SLX_barcode$SLX_ID, "_", patient_SLX_barcode$barcode.)
    patient_SLX_barcode <- patient_SLX_barcode[,"SLX_barcode"]
    print(paste0("Generated list of all patinets accross all cohorts that can be used for locus noise filter determination. Total: ", length(patient_SLX_barcode), " cases"))

    error_rate.raw <<- filter(combined.filtered, SLX_barcode %in% patient_SLX_barcode) %>% #case_or_control == "case"
      plyr::ddply(c("uniq_pos", "TRINUCLEOTIDE"), function(x){
        mut_sum <- sum(x$ALT_F) +sum(x$ALT_R)
        DP <- sum(x$DP)
        n_samples <- length(unique(x$SLX_barcode))

        n_samples_with_signal <- filter(x, ALT_F + ALT_R > 0)$SLX_barcode %>% unique() %>% length()

        data.frame(background_AF = mut_sum / DP, mut_sum, DP, n_samples, n_samples_with_signal)
      })

      save(error_rate.raw, file = "output_R/locus_error_rates.off_target.Rdata")
  }



  error_rate.raw$proportion <- error_rate.raw$n_samples_with_signal / error_rate.raw$n_samples

  error_rate.raw$LOCUS_NOISE.PASS <- error_rate.raw$proportion < proportion_of_controls & error_rate.raw$background_AF < max_background_mean_AF #

  nonzero_background_AFs <- filter(error_rate.raw, LOCUS_NOISE.PASS == TRUE, background_AF > 0)$background_AF

  non_zero_loci_pc <- signif((nrow(filter(error_rate.raw, background_AF > 0))/nrow(error_rate.raw))*100, digits = 3)
  locus_noise_fail_pc <- signif((nrow(filter(error_rate.raw, LOCUS_NOISE.PASS == F))/nrow(error_rate.raw))*100, digits = 3)

  if (on_target == TRUE){
    plot_dir <- "output_R/"
    print("running locally, will make a nice plot")
    error_rate.raw$COSMIC_plotting <- error_rate.raw$COSMIC
    error_rate.raw[error_rate.raw$COSMIC >20, "COSMIC_plotting"] <- "COSMIC > 20"
    p1 <- error_rate.raw %>%
      ggplot2::ggplot(aes(x = background_AF, fill = LOCUS_NOISE.PASS))+
      geom_histogram(bins = 100, position = "dodge")+
      scale_colour_discrete(name = "Locus noise pass")+
      #scale_y_log10()+
      scale_x_log10(breaks = c(1e-4, 1e-3, 1e-2, 1e-1))+
      #geom_vline(xintercept = locus_threshold, linetype = "dashed")+
      theme_bw()+
      labs(x = "Background plasma AF across all samples",
           y = "Frequency",
           title = paste(study, setting, "on-target, nonptspec data"),
           subtitle = paste0("Blacklisting of non-zero AF loci (representing ", non_zero_loci_pc, "% of data)\nBlacklisted loci (LOCUS_NOISE.FAIL) = " , locus_noise_fail_pc, "%\nsplit by COSMIC mutation status"))+
           facet_wrap(~COSMIC_plotting, scales = "free_y")

      error_rate.raw <- dplyr::select(error_rate.raw, -COSMIC_plotting)

    if (plot == TRUE){
    ggsave(p1, filename= paste0(plot_dir, "p1_", study, "_", setting, ".pdf"), width = 6, height = 4)
    }
  }


  combined.filtered$LOCUS_NOISE.PASS <- (combined.filtered$uniq_pos %in% filter(error_rate.raw, proportion < proportion_of_controls)$uniq_pos) &
                                              (combined.filtered$uniq_pos %in% filter(error_rate.raw, background_AF < max_background_mean_AF)$uniq_pos)

  return(combined.filtered)
}

annotate_with_both_strands <- function(curr){
  print("identifying mutations that were not actually supported by both F + R info now that we have size info")
  curr.mut <- filter(curr, mutant > 0,
                     BOTH_STRANDS == TRUE,
                     #PASS == TRUE,
                     LOCUS_NOISE.PASS == TRUE)

  n_molecules_in_both_strands <- plyr::ddply(curr.mut, "uniq_pos", function(x){
    data.frame(n_molecules = length(which(duplicated(x$size))))
  })

  print(paste("number of loci masquerading as being supported by both strands = ", nrow(filter(n_molecules_in_both_strands, n_molecules == 0))))

  # annotate the false F+R supported mutants as not passing the both_strands filter
  no_strands <- filter(n_molecules_in_both_strands, n_molecules == 0)$uniq_pos
  curr[curr$uniq_pos %in% no_strands, 'BOTH_STRANDS'] <- FALSE

  return(curr)
}

### TAPAS4 -------------------------------------
TAPAS4 <- function(data_dt.on_target,
                    cosmic_threshold = 0,
                    proportion_of_controls,
                    AF_threshold = 0.01, # don't use samples with >1% IMAF in their ptspec loci as controls at nonptspec sites
                    tumour_mutations,
                    error_file,
                    SLX_layout,
                    outlier_suppression_threshold,
                    max_DP,
                    max_background_mean_AF){

  print("running TAPAS4")
  #print(paste("cosmic mutations (data_dt.on_target):", table(data_dt.on_target$COSMIC_MUTATIONS)))
  starting_mutations <- length(unique(data_dt.on_target[data_dt.on_target$uniq_pos %in% tumour_mutations$uniq_pos,'uniq_pos']$uniq_pos))
  print(paste("number of mutations in data_dt.on_target that are also in tumour mutations", starting_mutations))

  # load input files
  combined <- parse(data_dt.on_target,
        tumour_mutations,
        error_file,
        SLX_layout,
        cosmic_threshold) %>%
        filter(DP < max_DP)

  #print(paste("cosmic mutations (combined, pre anotation with locus error rate):", table(combined$COSMIC_MUTATIONS)))
  print(paste("number of unique loci (post-parse)", length(unique(combined$uniq_pos))))

  # check error rate per locus
  print(paste("max background_AF allowed:", max_background_mean_AF))
  print(paste("max proportion of controls a locus can be positive in:", proportion_of_controls))
  combined <- annotate_with_locus_error_rate(combined.filtered = combined,
                                             proportion_of_controls,
                                             max_background_mean_AF,
                                             on_target = TRUE)

print(paste("number of unique loci (post-locus-noise)", length(unique(combined$uniq_pos))))

  if (grepl(study, pattern = "MELR")){
    print("MELR - not using dilution series for background determination as SNPs would be diluted in AF")
    nonptspec <- filter(combined, data == "nonptspec" & !grepl(";", Patient) & !grepl("dil", sample_name))
    ptspec <- filter(combined, data == "ptspec")

    combined.filtered <- rbind(nonptspec, ptspec)

    print("MELR - excluding patients who were not actually included on the panel")
    combined.filtered <- filter(combined.filtered, !grepl("MR1003", Patient), !grepl("MR1009", Patient))

    ## end of MELR-specific filters
  } else {
    combined.filtered <- combined
  }

  print(paste("using cosmic threshold from config of:", cosmic_threshold))
  combined.filtered <- filter(combined.filtered, data == "ptspec" | (data == "nonptspec" & COSMIC_MUTATIONS <= cosmic_threshold))

  # perform polish
  combined.polished <- combined.filtered %>%
                          both_strands() %>%
                          #background_polish(cluster_or_local = "local", p_value = outlier_suppression_threshold) %>%
                          as.data.table()

  ## annotate nonptspec with samples at risk of spillover
  combined.polished <- annotate_nonptspec_samples_to_skip(combined.polished, AF_threshold)

  print(paste("number of unique loci (post-polish)", length(unique(combined.polished$uniq_pos))))

  return(combined.polished)
}

## Parsing functions -------------------------------------
parse <- function(data_dt.on_target,
                  tumour_mutations,
                  error_file,
                  SLX_layout,
                  cosmic_threshold = 20,
                  skip_patient_from_background = NA){

  print(paste("setting cosmic threshold as", cosmic_threshold))
  data_dt.on_target$COSMIC <- data_dt.on_target$COSMIC_MUTATIONS > cosmic_threshold

  # update col_names of tumour mutations
  colnames(tumour_mutations)[grep(colnames(tumour_mutations), pattern = "^alt", ignore.case = T, perl = T)]  <- "ALT"
  colnames(tumour_mutations)[grep(colnames(tumour_mutations), pattern = "^ref", ignore.case = T, perl = T)]  <- "REF"

  columns_needed <- c("pt_uniq_pos", "tumour_AF", "Patient", "uniq_alt", "uniq_pos", "REF", "ALT")

  if (all(columns_needed %in% colnames(tumour_mutations)) == TRUE){
    print("tumour mutations columns correct")
  } else{
    print(paste("ERROR - Tumour mutations dataframe needs one of the following columns", columns_needed))
    print("Needs to be fixed by user.....")
    return(NULL)
  }

  colnames_SLX <- c("Study", "sample_name", "Patient", "sample_type", "input_into_library_ng",
           "SLX_barcode")

  if(all(colnames_SLX %in% colnames(SLX_layout)) == TRUE){
    print("SLX layout correct")

    data_dt.on_target <- as.data.frame(data_dt.on_target)

    # intersect SLX layout table with samples table
    data_dt.on_target <- annotate_with_SLX_table(data_dt.on_target, as.data.frame(SLX_layout),
                                                 colnames = c("Study", "sample_name", "Patient", "sample_type", "input_into_library_ng",
                                                              "SLX_barcode"))
  } else{
    print(paste("your SLX layout needs colnames:", colnames_SLX))
    return(NULL)
  }

  data_dt.on_target <- ignore_unexpected_variants(data_dt.on_target, tumour_mutations)

  # add mutation class
  tumour_mutations$mut_class <- paste0(tumour_mutations$REF, "/", tumour_mutations$ALT)
  print(paste("number of unique mutations in tumour_mutations file", length(unique(tumour_mutations$uniq_pos))))

  print("merging")
  ptspec <- semi_join(data_dt.on_target, tumour_mutations, by = "pt_uniq_pos")
  #print(paste0("Total rows of ptspec: ", nrow(ptspec)))
  print(paste("number of unique loci (ptspec) ", length(unique(ptspec$uniq_pos))))

  #saveRDS(data_dt.on_target, "on_target.annotated.RDS")

  nonptspec <- anti_join(data_dt.on_target, tumour_mutations, by = "pt_uniq_pos")
  print(paste("number of unique loci (nonptspec) ", length(unique(nonptspec$uniq_pos))))

  #print(paste0("Total rows of nonptspec: ", nrow(nonptspec)))

  ## annotate with tumour AF data
  print("annotating with tumour AF")
  ptspec <- inner_join(ptspec, tumour_mutations[!duplicated(tumour_mutations$pt_uniq_pos),
                                               c("pt_uniq_pos", "tumour_AF", "mut_class", "Patient")], by = "pt_uniq_pos")

  if (!is.na(skip_patient_from_background) == TRUE){
    print(paste0("skipping patients containing:", skip_patient_from_background))
    tumour_mutations <- tumour_mutations[!grepl(tumour_mutations$Patient, pattern = skip_patient_from_background), ]
  }

  print("generating nonptspec")
  nonptspec <- inner_join(nonptspec, tumour_mutations[!duplicated(tumour_mutations$uniq_pos),
                                                     c("uniq_pos", "tumour_AF", "mut_class", "Patient")], by = "uniq_pos")

  ptspec$data <- "ptspec"
  nonptspec$data <- "nonptspec"

  colnames(ptspec)[grepl(colnames(ptspec), pattern = "Patient")] <- c("Patient", "pt_mutation_belongs_to")
  colnames(nonptspec) <- colnames(ptspec)

  # n mutations before ptspec / nonptspec merge
  temp <- rbind(ptspec, nonptspec)
  print(paste("number of unique loci before combining ptspec + nonptspec", length(unique(temp$uniq_pos))))

  healthy_control_info <- filter(SLX_layout, case_or_control == "control")
  healthy_control_info <- healthy_control_info[,"SLX_barcode"]
  # ensure only loci that are interrogated in ptspec are looked for in the nonptspec
  print("ensure only loci that are interrogated in ptspec are looked for in the nonptspec + vice versa")
  ptspec <- ptspec[ptspec$uniq_pos %in% nonptspec$uniq_pos, ]
  nonptspec_cases <- nonptspec[!(nonptspec$SLX_barcode %in% healthy_control_info), ]
  print(paste0("Using following samples for nonptspec cases:", unique(nonptspec_cases$SLX_barcode)))
  nonptspec_cases <- nonptspec_cases[nonptspec_cases$uniq_pos %in% ptspec$uniq_pos, ]

  nonptspec_controls <- nonptspec[nonptspec$SLX_barcode %in% healthy_control_info, ]
  print(paste0("Using following samples for nonptspec controls:", unique(nonptspec_controls$SLX_barcode)))
  nonptspec_controls <- nonptspec_controls[nonptspec_controls$uniq_pos %in% ptspec$uniq_pos, ]

  combined <- rbind(ptspec, nonptspec_cases, nonptspec_controls)
  print(paste("number of unique loci after combining ptspec + nonptspec", length(unique(combined$uniq_pos))))

  print("combine classes")
  combined <- combine_classes_in_rds(combined)

  ## annotate with error rate
  SLX_layout <- data.table(SLX_layout)
  print("background error calculations")
  background_error <- calculate.background_error(error_file, SLX_layout)
  background_error <- add.missing_error_classes(background_error, trinucleotide_depth)
  background_error <- combine_classes(background_error)

  combined <- inner_join(combined, filter(background_error, case_or_control ==
                                            "case", data == "locus_noise.both_reads")[, c("TRINUCLEOTIDE", "mut_class", "background.mut_sum",
                                                                                         "background.DP", "background_AF")], by = c("TRINUCLEOTIDE", "mut_class"))

  print("Error rates successfully merged with locus specific data")

  combined$mut_sum <- combined$ALT_F + combined$ALT_R

  combined$sample_name <- paste0(combined$sample_name, " (", combined$pt_mutation_belongs_to, ")")
  print("Sample name column was added")
  #INVAR_input <- combined[, c("data", "Patient", "sample_name", "uniq_pos",
  #"DP", "mut_sum", "ALT_F", "ALT_R", "AF", "mut_class", "TRINUCLEOTIDE", "background_AF",
  #"tumour_AF", "COSMIC", "SLX_barcode", "uniq_alt", "uniq")]

  return(combined)
}


annotate_with_SLX_table <- function(data_dt.on_target, SLX_layout, colnames){
  data_dt.on_target <- left_join(data_dt.on_target, SLX_layout[,colnames], by = "SLX_barcode")

  # annotate on target file
  data_dt.on_target$uniq <- paste0(data_dt.on_target$uniq_pos, "_", data_dt.on_target$SLX_barcode)

  # add ptspec info
  data_dt.on_target$pt_uniq_pos <- paste0(data_dt.on_target$Patient, "_", data_dt.on_target$uniq_pos)
  #data_dt.on_target$pt_uniq_alt <- paste0(data_dt.on_target$pt_uniq_pos, "_", data_dt.on_target$REF, "/", data_dt.on_target$ALT)
  data_dt.on_target$uniq_alt <- paste0(data_dt.on_target$CHROM, ":", data_dt.on_target$POS, "_", data_dt.on_target$REF, "/", data_dt.on_target$ALT)

  return(data_dt.on_target)
}

## exclude mutations in patient-specific mpileup that are not observed in patient tumours
ignore_unexpected_variants <- function(merged, tumour_mutations){

  non_ref <- filter(merged, AF > 0)
  ref <- filter(merged, AF == 0)

  print(paste("non_ref", nrow(non_ref)))
  print(paste("ref", nrow(ref)))

  print('ignoring variants that do not correspond exactly with tumour variants')
  non_ref.removed <- non_ref[!(non_ref$uniq_alt %in% tumour_mutations$uniq_alt), ]
  non_ref <- non_ref[non_ref$uniq_alt %in% tumour_mutations$uniq_alt, ]

  print(paste("nrow of non_ref, removed", nrow(non_ref.removed)))
  print(paste("nrow of non_ref", nrow(non_ref)))

  #print('returning entirely patient-specific dataframe')
  merged <- rbind.fill(ref, non_ref) #chr3:51664915_C/G
  print("returning")
  print( paste("unique loci in merged after ignoring unexpecte variants:", length(unique(merged$uniq_pos))) )

  return(merged)
}

# combine complementary mutation classes
combine_classes_in_rds <- function(combined){
  print("combining complementary mutation classes and trinucleotides")
  combined.complement1 <- filter(combined, REF == "C" | REF== "T")
  combined.complement2 <- filter(combined, REF == "A" | REF== "G")

  if("TRINUCLEOTIDE" %in% colnames(combined.complement2)){
    print("swapping TRINUCLEOTIDE COLUMN, and reversing with reverse()")
    combined.complement2$TRINUCLEOTIDE <- chartr("ATGC","TACG", combined.complement2$TRINUCLEOTIDE)
    combined.complement2$TRINUCLEOTIDE <- reverse(as.character(combined.complement2$TRINUCLEOTIDE))
  } else{
    print("Error - your TRINUCLEOTIDE column will not be swapped. Check function code")
  }

  print("swapping mut_class column")
  combined.complement2$mut_class <- chartr("ATGC","TACG", combined.complement2$mut_class)

  combined1 <- rbind(combined.complement1, combined.complement2)

  return(combined1)
}

# identify samples with ctDNA >1% so that they are flagged
# later in the pipeline they are not used as patient-controls because of a
# potential risk of contamination between pt-spec and non-pt-spec loci
annotate_nonptspec_samples_to_skip <- function(combined.polished, AF_threshold = 0.01, plot_dir = "output_R/"){

  AF <<- filter(combined.polished, LOCUS_NOISE.PASS == TRUE, BOTH_STRANDS == TRUE) %>% #PASS == TRUE,
    plyr::ddply(c("sample_name", "data"), function(x){
      mut_sum <- sum(x$ALT_F) + sum(x$ALT_R)
      depth <- sum(x$DP)
      data.frame(mut_sum, depth, AF = mut_sum / depth)
    })

  ## any correlation between sample and it's nonptspec and ptspec AFs?
  AF$sample <- unlist(lapply(strsplit(AF$sample_name, split = " \\("), "[[", 1))
  AF$patient <- unlist(lapply(strsplit(AF$sample_name, split = " \\("), "[[", 2))

print(paste0("colnames of AF table :", colnames(AF)))

  test <- inner_join(filter(AF, data == "ptspec"),
                     filter(AF, data == "nonptspec")[,c(1,5,6)],
                     by = "sample")

  print(paste("filtering based on AF threshold", AF_threshold))
  low_samples <- filter(test, AF.x < AF_threshold)
  print(cor.test(low_samples$AF.x, low_samples$AF.y))

  do_not_use <- filter(AF, data == "ptspec", AF >AF_threshold)$sample

  AF$CONTAMINATION_RISK.PASS <- AF$sample %in% do_not_use
  combined.polished$CONTAMINATION_RISK.PASS <- !(unlist(lapply(strsplit(combined.polished$sample_name, split = " \\("), "[[", 1)) %in% do_not_use)

  return(combined.polished)
}

# For multi-panel data, check that the mutation being interrogated is on panel
ensure_on_panel <- function(curr, tumour_mutations){

  panel_key <- tumour_mutations[!duplicated(tumour_mutations$Patient),c("Patient", "panel_number")]

  curr <- left_join(curr, panel_key, by = "Patient")

  colnames(panel_key)[1] <- "pt_mutation_belongs_to"
  curr <- left_join(curr, panel_key, by = "pt_mutation_belongs_to")
  curr <- filter(curr, panel_number.x == panel_number.y)

  return(curr)
}

# Size filtering and annotation --------------------------------------
parse.size_ann <- function(size_ann, SLX_layout, tumour_mutations){
  size_ann$SLX_barcode <- gsub(size_ann$file_name, pattern = paste0(".*/", SLX_barcode_format), replacement = "\\1", perl = T)
  size_ann$uniq_pos <- gsub(size_ann$uniq, pattern = "(.*):[ACGT]:[ACGT]", perl = T, replacement = "\\1")
  colnames(size_ann)[2] <- "uniq_alt"

  size_ann$uniq <- paste0(size_ann$uniq_pos, "_", size_ann$SLX_barcode)
  size_ann$uniq_alt <- gsub(size_ann$uniq_alt, pattern = "(.*):([ACGT]):([ACGT])$", replacement = "\\1_\\2/\\3", perl = T)

  size_ann <- left_join(size_ann, combined[,c("uniq","data")] , by = "uniq")
  size_ann <- size_ann[!is.na(size_ann$data),]

  size_ann1 <- inner_join(size_ann, SLX_layout[,c("library_prep", "SLX_barcode")], by = "SLX_barcode")
  size_ann1$library_prep <- tolower(size_ann1$library_prep)
  return(size_ann1)
}

downsample_fragments <- function(size_ann, discrepant_uniq_pos, mpileup.DP, mpileup.mut_sum){
  print("starting downsample fragments function")
  temp <- size_ann[uniq_pos == discrepant_uniq_pos, ]

  temp.mutant <- filter(temp, mutant == 1)
  temp.wt <- filter(temp, mutant != 1)

  if(nrow(temp.mutant) != mpileup.mut_sum){
    print("discrepant mut_sum, going order the sizes by descending number of entries")
    temp.mutant.ordered <- temp.mutant[order(temp.mutant$size),]

    print(paste("nrow of temp.mutant = ", nrow(temp.mutant.ordered)))
    print("setting unpaired mutant reads not supported by mpileup to zero mutant status")
    temp.mutant.ordered[(mpileup.mut_sum + 1) : nrow(temp.mutant.ordered), 'mutant'] <- 0

    temp.mutant <- temp.mutant.ordered
  }

  temp <- rbind(temp.wt, temp.mutant)

  temp.mutant <- filter(temp, mutant == 1)
  temp.wt <- filter(temp, mutant != 1)

  # in some cases, the number of reads from the size CSV differs from the number of mpileup reads
  if(nrow(temp) > mpileup.DP){
    print("python size data has more rows than mpileup data")
    print("discrepant DP, downsampling wildtypes from the mpileup data")
    wt_reads <- mpileup.DP - mpileup.mut_sum

    print(paste("wt reads:", wt_reads))
    print(paste("nrows of temp wt: ", nrow(temp.wt)))

    temp.wt <- temp.wt[sample(nrow(temp.wt), wt_reads, replace = F),]
  } else if (nrow(temp) < mpileup.DP){
    print("python size data has FEWER rows than mpileup data")
    print("discrepant DP, will resample (with replacement) aka stretch mpileup data for wt data (inconsequential as wt reads are not used later in pipeline)")
    wt_reads <- mpileup.DP - mpileup.mut_sum

    print(paste("wt reads:", wt_reads))
    print(paste("nrows of temp wt: ", nrow(temp.wt)))

    temp.wt <- temp.wt[sample(nrow(temp.wt), wt_reads, replace = T),]
  }

  print("rbinding temp.wt and temp.mutant")

  output <- rbind(temp.wt, temp.mutant)
  return(output)
}

equalise_size_counts <- function(size_ann.csv = "/Users/wan01/groups/Research/nrlab/group_folders/TAPAS/input/MELR/size.split/SLX-11848_D701_D501.realign.recal.BL.connor.f0.9_s2.bam.inserts_for_annotation.csv",combined.polished){
  print("reading size_ann")
  curr_SLX_barcode <- gsub(size_ann.csv, pattern = SLX_barcode_format, perl = T, replacement = "\\1")
  curr_SLX_barcode <- unlist(lapply(strsplit(curr_SLX_barcode, split = "\\."),"[[",1))
  size_ann <- read.csv(file = size_ann.csv, col.names = c("ALT", "uniq", "mutant", "size"))

  print(paste("curr SLX barcode: ", curr_SLX_barcode))
  print(unique(combined.polished$SLX_barcode))

  if ((curr_SLX_barcode %in% unique(combined.polished$SLX_barcode)) == FALSE){
    print("This SLX barcode was not observed in your mpileup (combined) df. size annotation cannot be performed")
    return()
  }

  curr <- filter(combined.polished, SLX_barcode == curr_SLX_barcode)

  ## fix DP count in curr
  curr$DP <- curr$ALT_F + curr$ALT_R + curr$REF_F + curr$REF_R

  size_ann <- data.table(size_ann)
  size_ann$mutant[size_ann$mutant == "."] <- 0
  size_ann$mutant <- as.numeric(as.character(size_ann$mutant))

  size_ann$SLX_barcode <- curr_SLX_barcode
  size_ann$mut_class <- gsub(size_ann$uniq, pattern = ".*:([ACGT]):([ACGT])$", replacement = "\\1/\\2",
                             perl = T)
  size_ann$uniq_pos <- gsub(size_ann$uniq, pattern = "(.*):([ACGT]):([ACGT])$", replacement = "\\1",
                            perl = T)

  size_ann$REF <- gsub(size_ann$uniq, pattern = ".*:([ACGT]):([ACGT])$", replacement = "\\1", perl = T)

  ## reduce 12 mut classes to 6
  # columns of this data.frame = base, chr:pos:ref:alt, mutant (boolean), size (bp)
  size_ann <- combine_classes_in_rds(size_ann) %>% data.table()

  print("collapse size ann table")
  size_ann.collapsed <- size_ann[, .(mut_sum.size =sum(mutant), DP.size=.N) , by=.(uniq, uniq_pos, mut_class, size)]

  size_ann.summary <- size_ann.collapsed[, .(DP.size = sum(DP.size),
                                              mut_sum.size = sum(mut_sum.size)),
                                         by = list(uniq_pos, mut_class)]

  test <- left_join(curr, size_ann.summary, by = c("uniq_pos", "mut_class"))

  discrepant_DP <- filter(test, DP != DP.size |  ALT_F + ALT_R != mut_sum.size)
  correct <- filter(test, DP == DP.size & ALT_F + ALT_R == mut_sum.size)

  if (length(which(is.na(test$mut_sum.size))) > 0){
    print("stop; the size test df has loci which do not have size data")
    return(NULL)
  }

  ## correct the raw size df - only AF == 0
  discrepant_DP.zero <- filter(discrepant_DP, AF == 0) %>% data.table()
  discrepant_DP.nonzero <- filter(discrepant_DP, AF > 0)

  size_ann.discrepant.zero.fixed <- plyr::ddply(discrepant_DP.zero, "uniq_pos", function(x){
    print(paste("running...", unique(x$uniq_pos)))

    print("AF = 0 , just taking the top n molecule sizes at this locus")
    size_ann[uniq_pos == unique(x$uniq_pos)][1:unique(x$DP),]

  })

  ## set size_ann.discrepant.zero.fixed mutant counts to zero. This is because the mutant reads identified from size info
  # are not subject to BQ and MQ filters
  size_ann.discrepant.zero.fixed$mutant <- 0

  ## correct loci with AF > 0
  size_ann.discrepant.nonzero.fixed <- plyr::ddply(discrepant_DP.nonzero, "uniq_pos", function(x){
    print(paste("running...", unique(x$uniq_pos)))
      print("AF > 0, going to downsample loci in the WT")

      discrepant_uniq_pos <- unique(x$uniq_pos)
      mpileup.DP <- x$DP
      mpileup.mut_sum <- sum(x$ALT_F) + sum(x$ALT_R)

      print(paste("mpileup.DP = ", mpileup.DP))
      print(paste("mpileup.mut_sum = ", mpileup.mut_sum))

      downsample_fragments(size_ann = size_ann, discrepant_uniq_pos, mpileup.DP, mpileup.mut_sum)
  })

  print("1st rbind")
  size_ann.discrepant.fixed <- rbind(size_ann.discrepant.zero.fixed, size_ann.discrepant.nonzero.fixed)

  print("2nd rbind")
  size_ann.fixed <- rbind(size_ann.discrepant.fixed, size_ann[size_ann$uniq_pos %in% correct$uniq_pos,])

  print("merging corrected size df and original combined df")
  size_ann.fixed <- size_ann.fixed %>% dplyr::select(-SLX_barcode)
  combined.size_ann <- inner_join(curr, size_ann.fixed[,3:6], by = c("uniq_pos", "mut_class"))

  return(combined.size_ann)
}


# insert size code for visualisation of fragment sizes
load.insert_sizes <- function(file){
  print("reading in insert_sizes file, please wait")
  insert_sizes <- readr::read_delim(file = file, delim = "\t", col_names = F)#
  colnames(insert_sizes) <- c("name", "count")

  print("annotating with additional columns")
  # split out relevant info from the name
  insert_sizes$length <- gsub(insert_sizes$name, perl = T, pattern = ".*:([0-9]+)$", replacement = "\\1")
  insert_sizes$uniq_pos <- gsub(insert_sizes$name, perl = T, pattern = ".*_(chr.+:[0-9]+):.*", replacement = "\\1")
  insert_sizes$mut <- !grepl(insert_sizes$name, pattern = "ref")

  print("obtaining SLX and barcode information")
  insert_sizes$SLX <- gsub(insert_sizes$name, perl = T, pattern = ".*(SLX-[0-9]+)_.*", replacement = "\\1")
  insert_sizes$BARCODE <- gsub(insert_sizes$name, perl = T, pattern = ".*(D[0-9]+_D[0-9]+).*", replacement = "\\1")
  insert_sizes$SLX_barcode <- paste0(insert_sizes$SLX, "_", insert_sizes$BARCODE)

  insert_sizes$length <- as.numeric(insert_sizes$length)
  insert_sizes$uniq <- paste0(insert_sizes$uniq_pos, "_", insert_sizes$SLX_barcode)
  return(insert_sizes)
}
get.size <- function(insert_sizes.ptspec, round){
  #print("plotting size distribution and enrichment ratio")
  insert_sizes.ptspec$length.round <- round_any(insert_sizes.ptspec$length, round)

  # total number of read families
  size.totals <- ddply(insert_sizes.ptspec, c('mut', "sample_type", "case_or_control", "Study", "data"), function(x){
    data.frame(total = nrow(x))
  })

  size.sum <- ddply(insert_sizes.ptspec, c('mut', "sample_type", "case_or_control", "Study", 'length.round', "data"), function(x){
    data.frame(mut_sum = nrow(x))
  })

  size <- merge(size.totals, size.sum, by = c('mut', "sample_type", "case_or_control", "Study", "data"))
  names(size)[names(size) == 'mut_sum'] <- 'count'
  size$proportion <- as.numeric(size$count)/size$total

  return(size)
}

get.size.characterisation <- function(path_to_inserts_txt, combined.polished, SLX_layout, rounding_for_plot = 5){
  insert_sizes <- load.insert_sizes(file = path_to_inserts_txt)

  # annotate with size
  insert_sizes.stretch <- insert_sizes[rep(seq(1, nrow(insert_sizes)), insert_sizes$count), ]
  insert_sizes.stretch <- data.table(insert_sizes.stretch)

  insert_sizes.summary <- insert_sizes.stretch[, .(median = median(length),
                                                   max = max(length),
                                                   min = min(length),
                                                   sum = length(length)),
                                               by = list(SLX_barcode, uniq_pos, mut)]
  print("Fragment analysis of ptspec")
  loci <- (filter(combined.polished, data == "ptspec"))$uniq

  insert_sizes.ptspec <- insert_sizes.stretch[insert_sizes.stretch$uniq %in% loci, ]
  insert_sizes.ptspec <- left_join(insert_sizes.ptspec, SLX_layout[ , c("SLX_barcode", "sample_type", "Study", "case_or_control")], by = "SLX_barcode")
  insert_sizes.ptspec$data <- "ptspec"

  print("Fragment analysis of nonptspec")
  nonptspec <- filter(combined.polished, data == "nonptspec")
  insert_sizes.nonptspec <- insert_sizes.stretch[insert_sizes.stretch$uniq %in%
                                                   nonptspec$uniq, ]

  insert_sizes.nonptspec <- left_join(insert_sizes.nonptspec, SLX_layout[ , c("SLX_barcode", "sample_type", "Study", "case_or_control")], by = "SLX_barcode")
  insert_sizes.nonptspec$data <- "nonptspec"

  print(paste("binning the data into bins of size", rounding_for_plot))
  size.ptspec <- get.size(insert_sizes.ptspec, round = rounding_for_plot)
  size.nonptspec <- get.size(insert_sizes.nonptspec, round = rounding_for_plot)
  print("Merging ptspec and nonptspec dataframes")
  size.combined <- rbind(size.ptspec, size.nonptspec)
  size.combined$mut <- factor(size.combined$mut, levels = c(TRUE, FALSE))

  print(ggplot(size.combined, aes(x = length.round, y = proportion, fill = mut)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Fragment size in bp",
         y = "Proportion",
         title = "Fragment comparison",
         subtitle = paste(study, setting)) +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14,
                                    face = "bold"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12,
                                      face = "bold")) +
    scale_x_continuous(limits = c(0, 400))+
    facet_wrap(data~case_or_control, scales = "free_y"))
  ggsave(paste0(plot_dir, "p4_", study, "_", setting, ".pdf"), width = 8, height = 7)

  print("returning an un-binned dataframe of sizes and mutation status")
  size.ptspec <- get.size(insert_sizes.ptspec, round = 1)
  size.nonptspec <- get.size(insert_sizes.nonptspec, round = 1)
  size.combined <- rbind(size.ptspec, size.nonptspec)
  size.combined$mut <- factor(size.combined$mut, levels = c(TRUE, FALSE))

  return(size.combined)
}

## the size dataframe is an aggregate of all the samples, so
## when assessing the detection of a particular sample, exclude itself
## from the aggregated df so that it is not circular
leave_one_out.size <- function(raw.filtered, size_characterisation){
  print("ddply")
  mutant_reads_to_subtract <-
    filter(raw.filtered, mutant == TRUE) %>% ddply("size", function(x){
    data.frame(to_subtract = nrow(x)) })

  if (nrow(mutant_reads_to_subtract) > 0){
    print("mutant reads present")
    #print("ptspec only")
    #size_characterisation <- filter(size_characterisation, data == "ptspec")

    print("joining")
    size_characterisation <- left_join(size_characterisation, mutant_reads_to_subtract, by = "size")
    size_characterisation[is.na(size_characterisation$to_subtract),"to_subtract"] <- 0

    ## do not subtract mutant reads from the wild-type bins
    size_characterisation[size_characterisation$mut == FALSE, 'to_subtract'] <- 0

    ## update columns after subtraction
    print("leave_one_out - updating columns 1/2")
    size_characterisation$count <- pmax(size_characterisation$count - size_characterisation$to_subtract, 0)
    size_characterisation$total <- pmax(size_characterisation$total - size_characterisation$to_subtract, 0)
    size_characterisation$proportion <- size_characterisation$count/size_characterisation$total

    print("leave_one_out - updating columns 2/2")
    total_mut_reads <- sum(filter(size_characterisation, mut == T)$count)
    size_characterisation[size_characterisation$mut == T, 'total'] <- total_mut_reads

    size_characterisation.leave_one_out <- select(size_characterisation, -to_subtract)

    print("finished leave one out")
    return(size_characterisation.leave_one_out)
  } else{
    print("no mutant reads, not doing leave one out subtraction")
    return(size_characterisation)
  }
}


# Error rate functions --------------------------------------
# calculate background error - used in pipeline
calculate.background_error <- function(error_file, SLX_layout, exclude_PPC = F){  # PPC = control sample naming in NR lab

  print("loading error Rdata file")
  load(file = error_file, .GlobalEnv)
  print(paste("names of objects in error rate list: ", names(error_rates)))

  error_rate.one_read <- error_rates[[1]]
  error_rate.locus_noise <- error_rates[[2]]
  #error_rate.polished_only <- error_rates[[3]]
  error_rate.both_reads <- error_rates[[3]]
  error_rate.locus_noise.both_reads <- error_rates[[4]]

  error_rate.one_read$data <- "one_read"
  error_rate.locus_noise$data <- "locus_noise"
  error_rate.both_reads$data <- "both_reads"
  #error_rate.polished_only$data <- "polished_only"
  error_rate.locus_noise.both_reads$data <- "locus_noise.both_reads"

  # If the locus noise dataframe is completely clean - this means you have either a) set the locus noise threshold too low or b) not run enough control samples
  if(sum(error_rate.locus_noise$mut_sum) == 0){
    print("LOCUS NOISE mutant sum is zero")
    print("this is likely due to insufficient controls being run")
    print("either run more controls or override this with a higher locus noise filter (beware higher background rates)")

    return(NULL)
  }

  error_rates.both <- rbind(error_rate.one_read,
                            error_rate.locus_noise,
                            error_rate.both_reads,
                            #error_rate.polished_only,
                            error_rate.locus_noise.both_reads)

  error_rates.both$SLX_barcode <- paste0(error_rates.both$SLX, "_", error_rates.both$BARCODE)

  error_rates.both <- left_join(error_rates.both, SLX_layout[,c("sample_name", "case_or_control", "SLX_barcode")], by = "SLX_barcode")

  if( exclude_PPC == TRUE){
    print("excluding PPC")
    error_rates.both <- filter(error_rates.both, !grepl("PPC", sample_name))
  }

  # Calculate  error rate
  print("calculating background error rates")
  background_error <- plyr::ddply(error_rates.both, c("REF", "ALT", "TRINUCLEOTIDE", "case_or_control", "data"), function(x){
    mut_sum=sum(x$mut_sum)
    total_DP=sum(x$total_DP)
    background_AF= mut_sum / total_DP
    data.frame(mut_sum, total_DP, background_AF)
  })

  trinucleotide_depth <<- plyr::ddply(background_error, c("TRINUCLEOTIDE", "case_or_control", "data"), function(x){
    DP = sum(x$total_DP)
    TRINUCLEOTIDE = unique(x$TRINUCLEOTIDE)
    data.frame(TRINUCLEOTIDE, DP)
  })

  background_error2 <- plyr::ddply(background_error, c("TRINUCLEOTIDE", "case_or_control", "data"), function(x){
    DP = sum(x$total_DP)

    plyr::ddply(filter(x, ALT != "."), c("TRINUCLEOTIDE", "case_or_control", "data", "REF", "ALT"), function(x){
      mut_sum_total = x$mut_sum
      background_AF = (x$mut_sum) / DP
      data.frame(mut_sum_total, DP, background_AF)
    })
  })

  # count the number of missing classes
  print(
    plyr::ddply(background_error2, c("case_or_control", "data"), function(x){
      data.frame(contexts_represented = length(unique(paste(x$TRINUCLEOTIDE, "_", x$ALT))))
    }))

  # workaround to avoid mut_sum breaking the function :(
  if("mut_sum" %in% colnames(background_error2) == TRUE){
    background_error2 <- dplyr::select(background_error2, -mut_sum, -total_DP)
  }

  return(background_error2)
}

# other background error calculation functions
calculate.background_error_with_downsampling <- function(error_file){

  print("loading error Rdata file ")
  load(file = error_file, .GlobalEnv)
  error_rate.one_strand <- error_rates[[1]]
  error_rate.both_strands <- error_rates[[2]]

  error_rate.one_strand$data <- "one_strand"
  error_rate.both_strands$data <- "both_strands"

  error_rates.both <- rbind(error_rate.one_strand, error_rate.both_strands)

  error_rates.both$SLX_barcode <- paste0(error_rates.both$SLX, "_", error_rates.both$BARCODE)

  error_rates.both <- left_join(error_rates.both, SLX_layout[,c("sample_name", "case_or_control", "SLX_barcode")], by = "SLX_barcode")

  #error_rates.both <- error_rates.both[!grepl(error_rates.both$sample_name, pattern =  "Control"),]

  case <- filter(error_rates.both, case_or_control == "case")
  control <- filter(error_rates.both, case_or_control == "control_negative")

  case.sampled <- case[sample(nrow(case), nrow(control), replace = T),]

  error_rates.both.sampled <- rbind(case.sampled, control)

  # Calculate  error rate
  print("calculating background error rates")
  background_error <- plyr::ddply(error_rates.both.sampled, c("REF", "ALT", "TRINUCLEOTIDE", "case_or_control", "data"), function(x){
    mut_sum=sum(x$mut_sum)
    total_DP=sum(x$total_DP)
    background_AF= mut_sum / total_DP
    data.frame(mut_sum, total_DP, background_AF)
  })

  trinucleotide_depth <<- plyr::ddply(background_error, c("TRINUCLEOTIDE", "case_or_control", "data"), function(x){
    DP = sum(x$total_DP)
    TRINUCLEOTIDE = unique(x$TRINUCLEOTIDE)
    data.frame(TRINUCLEOTIDE, DP)
  })

  background_error2 <- plyr::ddply(background_error, c("TRINUCLEOTIDE", "case_or_control", "data"), function(x){
    DP = sum(x$total_DP)

    plyr::ddply(filter(x, ALT != "."), c("TRINUCLEOTIDE", "case_or_control", "data", "REF", "ALT"), function(x){
      mut_sum_total = x$mut_sum
      background_AF = (x$mut_sum) / DP
      data.frame(mut_sum_total, DP, background_AF)
    })
  })

  # count the number of missing classes
  print(
    plyr::ddply(background_error2, c("case_or_control", "data"), function(x){
      data.frame(contexts_represented = length(unique(paste(x$TRINUCLEOTIDE, "_", x$ALT))))
    }))

  #background_error2 <- select(background_error2, - mut_sum, - total_DP)
  return(background_error2[,1:8])
}

calculate.background_error_combined_files <- function(error_rates){

  error_rate.one_strand <- error_rates[[1]]
  error_rate.both_strands <- error_rates[[2]]

  error_rate.one_strand$data <- "one_strand"
  error_rate.both_strands$data <- "both_strands"

  error_rates.both <- rbind(error_rate.one_strand, error_rate.both_strands)

  error_rates.both$SLX_barcode <- paste0(error_rates.both$SLX, "_", error_rates.both$BARCODE)

  error_rates.both <- left_join(error_rates.both, SLX_layout[,c("sample_name", "case_or_control", "SLX_barcode")], by = "SLX_barcode")

  # Calculate background error rate
  background_error <- plyr::ddply(error_rates.both, c("REF", "ALT", "TRINUCLEOTIDE", "case_or_control", "data"), function(x){
    mut_sum=sum(x$mut_sum)
    total_DP=sum(x$total_DP)
    background_AF= mut_sum / total_DP
    data.frame(mut_sum, total_DP, background_AF)
  })

  trinucleotide_depth <<- plyr::ddply(background_error, c("TRINUCLEOTIDE", "case_or_control", "data"), function(x){
    DP = sum(x$total_DP)
    TRINUCLEOTIDE = unique(x$TRINUCLEOTIDE)
    data.frame(TRINUCLEOTIDE, DP)
  })

  background_error2 <- plyr::ddply(background_error, c("TRINUCLEOTIDE", "case_or_control", "data"), function(x){
    DP = sum(x$total_DP)

    plyr::ddply(filter(x, ALT != "."), c("TRINUCLEOTIDE", "case_or_control", "data", "REF", "ALT"), function(x){
      mut_sum_total = x$mut_sum
      background_AF = (x$mut_sum) / DP
      data.frame(mut_sum_total, DP, background_AF)
    })
  })

  # count the number of missing classes
  print(
    plyr::ddply(background_error2, c("case_or_control", "data"), function(x){
      data.frame(contexts_represented = length(unique(paste(x$TRINUCLEOTIDE, "_", x$ALT))))
    }))

  return(background_error2)
}
calculate.background_error_tag_plasma_seq <- function(error_file){

  print("loading error Rdata file ")
  load(file = error_file, .GlobalEnv)
  error_rate.one_strand <- error_rates[[1]]
  error_rate.both_strands <- error_rates[[2]]

  error_rate.one_strand$data <- "one_strand"
  error_rate.both_strands$data <- "both_strands"

  error_rates.both <- rbind(error_rate.one_strand, error_rate.both_strands)

  error_rates.both$SLX_barcode <- paste0(error_rates.both$SLX, "_", error_rates.both$BARCODE)

  error_rates.both <- left_join(error_rates.both, SLX_layout[,c("sample_name", "case_or_control", "SLX_barcode", "library_prep")], by = "SLX_barcode")

  # Calculate  d error rate
  print("calculating background error rates")
  background_error <- plyr::ddply(error_rates.both, c("REF", "ALT", "TRINUCLEOTIDE", "data", "library_prep"), function(x){
    mut_sum=sum(x$mut_sum)
    total_DP=sum(x$total_DP)
    background_AF= mut_sum / total_DP
    data.frame(mut_sum, total_DP, background_AF)
  })

  trinucleotide_depth <<- plyr::ddply(background_error, c("TRINUCLEOTIDE", "library_prep", "data"), function(x){
    DP = sum(x$total_DP)
    TRINUCLEOTIDE = unique(x$TRINUCLEOTIDE)
    data.frame(TRINUCLEOTIDE, DP)
  })

  background_error2 <- plyr::ddply(background_error, c("TRINUCLEOTIDE", "library_prep", "data"), function(x){
    DP = sum(x$total_DP)

    plyr::ddply(filter(x, ALT != "."), c("TRINUCLEOTIDE", "library_prep", "data", "REF", "ALT"), function(x){
      mut_sum_total = x$mut_sum
      background_AF = (x$mut_sum) / DP
      data.frame(mut_sum_total, DP, background_AF)
    })
  })

  # count the number of missing classes
  print(
    plyr::ddply(background_error2, c("library_prep", "data"), function(x){
      data.frame(contexts_represented = length(unique(paste(x$TRINUCLEOTIDE, "_", x$ALT))))
    }))

  #background_error2 <- select(background_error2, - mut_sum, - total_DP)
  return(background_error2)
}

# add missing error classes
add.missing_error_classes <- function(background_error, trinucleotide_depth){

  # adding in missing mutation classes
  background_error$mut_class <- paste0(background_error$REF, "/", background_error$ALT)
  background_error$class_context <- paste0(background_error$TRINUCLEOTIDE, "_", background_error$mut_class)

  unique_classes <- unique(background_error$class_context)
  combined.output <- data.frame()

  missing_classes <- plyr::ddply(background_error, c("case_or_control", "data"), function(x){
    c_o_c <- unique(x$case_or_control)
    strandiness <- unique(x$data)
    for (y in unique_classes){
      if (nrow(filter(x, class_context == y)) == 0){

        print('adding row because class was not observed')
        context <- lapply(strsplit(y, "_"), "[[", 1)
        mut_class <- lapply(strsplit(y, "_"), "[[", 2)
        ref <- lapply(strsplit(as.character(mut_class), "/"), "[[", 1)
        alt <- lapply(strsplit(as.character(mut_class), "/"), "[[", 2)
        temp <- data.frame(context, c_o_c, strandiness, alt, ref, 0, mut_class, y)
        colnames(temp) <- c("TRINUCLEOTIDE", "case_or_control", "data", "ALT", "REF", "mut_sum_total", "mut_class","class_context")
        combined.output <- rbind(temp, combined.output)

      } else{
        print('class observed, skip')
      }
    }
    missing_classes <- combined.output
  })

  if (nrow(missing_classes) == 0){
    print("All classes are observed")
  } else {
    missing_classes <- left_join(missing_classes, trinucleotide_depth, by = c("case_or_control", "data", "TRINUCLEOTIDE"))
    missing_classes$background_AF <- missing_classes$mut_sum_total / missing_classes$DP
  }

  # now combine background_error and missing classes to final dataframe
  background_error_final <- rbind(background_error, missing_classes)

  return(background_error_final)
}

# combine complementary classes and trinucleotide contexts i.e. A>T etc.
combine_classes <- function(background_error){
  background_error.complement1 <- filter(background_error, REF == "C" | REF== "T")
  background_error.complement2 <- filter(background_error, REF == "A" | REF== "G")

  background_error.complement2$TRINUCLEOTIDE <- chartr("ATGC","TACG", background_error.complement2$TRINUCLEOTIDE)

  # flip context (and include a reverse because of the opposite strand)
  background_error.complement2$TRINUCLEOTIDE <- reverse(as.character(background_error.complement2$TRINUCLEOTIDE))
  background_error.complement2$REF <- chartr("ATGC","TACG", background_error.complement2$REF)
  background_error.complement2$ALT <- chartr("ATGC","TACG", background_error.complement2$ALT)
  background_error.complement2$mut_class <- chartr("ATGC","TACG", background_error.complement2$mut_class)
  background_error.complement2$class_context <- paste0(background_error.complement2$TRINUCLEOTIDE, "_", background_error.complement2$mut_class)

  background_error.combined <- rbind(background_error.complement1, background_error.complement2)

  background_error.combined.average <- plyr::ddply(background_error.combined, c("TRINUCLEOTIDE", "case_or_control", "data", "REF", "ALT", "mut_class"), function(x){
    mut_sum_total = sum(x$mut_sum_total)
    DP = sum(x$DP)

    data.frame(mut_sum_total, DP, background_AF = mut_sum_total/DP)
  })

  print("updating colnames - adding 'background' into mut_sum and DP colnames")
  colnames(background_error.combined.average)[c(7,8)] <- c("background.mut_sum", "background.DP")
  print("Columns have been updated")

  return(background_error.combined.average)
}


# INVAR/GLRT detection functions ------------------------
calculate_IMAFv2 <- function(raw, is.blood_spot = FALSE){
  # raw is a loaded combined.polished.size.os.Rdata file, same as what is being used by GLRT
  # is.blood_spot is defined in config2.R when low input samples are used (omits the outlier suppression step)
  print("calculating IMAFv2 i.e. BS_uAF")
  raw_flat <- select(raw, -size, -mutant) %>% distinct()

  if (is.blood_spot == TRUE){
    print("blood spot data, no need to use outlier suppression when depth is either 1 or 0")
    raw_flat <- filter(raw_flat,
                     LOCUS_NOISE.PASS == TRUE,
                     BOTH_STRANDS == TRUE)
  } else {
    raw_flat <- filter(raw_flat,
                     LOCUS_NOISE.PASS == TRUE,
                     BOTH_STRANDS == TRUE,
                     PASS == TRUE)
  }
  ##get mutation rate per trinucleotide context
  curr_summary <- plyr::ddply(raw_flat, c("mut_class", "TRINUCLEOTIDE"), function(x){
    total_DP <- sum(x$REF_F) + sum(x$REF_R) + sum(x$ALT_F) + sum(x$ALT_R)
    mut_sum <- sum(x$mut_sum)
    mean_AF <- weighted.mean(x$AF, x$DP)
    background_AF.trinucleotide <- unique(x$background_AF)

    data.frame(total_DP, mut_sum, mean_AF, background_AF.trinucleotide,
               mean_AF.BS_trinucleotide = mean_AF - background_AF.trinucleotide)
  })

  curr_summary$mean_AF.BS_trinucleotide <- pmax(curr_summary$mean_AF.BS_trinucleotide, 0)

  IMAFv2 <- signif(weighted.mean(curr_summary$mean_AF.BS_trinucleotide, curr_summary$total_DP), 4)
  print(paste("IMAFv2 =", IMAFv2))
  return(IMAFv2)
}

get.INVAR_score <- function(INVAR_score.path, SLX_layout,
                            adjust = T,
                            outlier_suppression = 0.05,
                            filter.pattern = NA,
                            CONTAMINATION_RISK_THRESHOLD = 0.01,
                            INVAR_score_specificity = 0.95){

  if(is.null(nrow(SLX_layout))){
    print("your SLX_layout is not a dataframe, please provide a dataframe (and not a string)")
    return(NULL)
  }

  print(paste0("Adjusting INVAR scores = ", adjust, ",outlier suppression = ", outlier_suppression))
  #print(paste("recalculate contamination risk?", !skip_contamination_reassessment))
  INVAR_scores <- read.table(file = INVAR_score.path, col.names = c("FILE_NAME", "OUTPUT"))

  # handle mixture experiment naming
  if (all(grep(INVAR_scores$OUTPUT, pattern = "; MR"))== TRUE){
    print("mixture samples detected based on MELR naming")
    INVAR_scores$OUTPUT <- gsub(INVAR_scores$OUTPUT, pattern = "; MR", replacement = "_MR")
  }

  paste("parsing col names")
  INVAR_scores$GLRT_setting <- unlist(lapply(strsplit(as.character(INVAR_scores$OUTPUT), split = ","), "[[",1))

  INVAR_scores$OUTPUT <- unlist(lapply(strsplit(as.character(INVAR_scores$OUTPUT), split = ": "), "[[",2))

  INVAR_scores$INVAR_SCORE <- unlist(lapply(strsplit(as.character(INVAR_scores$OUTPUT), split = "\\ "), "[[",1)) %>% as.numeric()
  INVAR_scores$DP <- unlist(lapply(strsplit(as.character(INVAR_scores$OUTPUT), split = "\\ "), "[[",2)) %>% as.numeric()
  INVAR_scores$AF_p <- unlist(lapply(strsplit(as.character(INVAR_scores$OUTPUT), split = "\\ "), "[[",3)) %>% as.numeric()
  INVAR_scores$IMAF <- unlist(lapply(strsplit(as.character(INVAR_scores$OUTPUT), split = "\\ "), "[[",4)) %>% as.numeric()
  INVAR_scores$patient <- unlist(lapply(strsplit(as.character(INVAR_scores$OUTPUT), split = "\\ "), "[[",5))
  INVAR_scores$pt_mut_belongs_to <- unlist(lapply(strsplit(as.character(INVAR_scores$OUTPUT), split = "\\ "), "[[",6))
  INVAR_scores$POLISH.PASS <- unlist(lapply(strsplit(as.character(INVAR_scores$OUTPUT), split = "\\ "), "[[",7))
  INVAR_scores$BOTH_STRANDS.PASS <- unlist(lapply(strsplit(as.character(INVAR_scores$OUTPUT), split = "\\ "), "[[",8))
  INVAR_scores$LOCUS_NOISE.PASS <- unlist(lapply(strsplit(as.character(INVAR_scores$OUTPUT), split = "\\ "), "[[",9))

  if (length(strsplit(as.character(INVAR_scores$OUTPUT), split = "\\ ")[[1]]) >= 15){
    print("Keep on splitting into further columns")
    INVAR_scores$smooth <- (unlist(lapply(strsplit(as.character(INVAR_scores$OUTPUT), split = "\\ "), "[[",10)))
    INVAR_scores$called_against <- unlist(lapply(strsplit(as.character(INVAR_scores$OUTPUT), split = "\\ "), "[[",11))
    INVAR_scores$CONTAMINATION_RISK.PASS <- (unlist(lapply(strsplit(as.character(INVAR_scores$OUTPUT), split = "\\ "), "[[",12)))
    INVAR_scores$outlier_suppression <- unlist(lapply(strsplit(as.character(INVAR_scores$OUTPUT), split = "\\ "), "[[",13))
    INVAR_scores$null_likelihood <- as.numeric(unlist(lapply(strsplit(as.character(INVAR_scores$OUTPUT), split = "\\ "), "[[",14)))
    INVAR_scores$alt_likelihood <- unlist(lapply(strsplit(as.character(INVAR_scores$OUTPUT), split = "\\ "), "[[",15))
    INVAR_scores$mutant_present <- unlist(lapply(strsplit(as.character(INVAR_scores$OUTPUT), split = "\\ "), "[[",16))

    print("filtering for one outlier_suppression condition")
    print(nrow(INVAR_scores))
    print(unique(INVAR_scores$outlier_suppression))
    INVAR_scores <- INVAR_scores[grepl(INVAR_scores$outlier_suppression, pattern = outlier_suppression),]
    print(nrow(INVAR_scores))

    INVAR_scores$n_mutations <- unlist(lapply(strsplit(as.character(INVAR_scores$OUTPUT), split = "\\ "), "[[",17))
    INVAR_scores$mut_sum <- unlist(lapply(strsplit(as.character(INVAR_scores$OUTPUT), split = "\\ "), "[[",18))
    INVAR_scores$iteration <- unlist(lapply(strsplit(as.character(INVAR_scores$OUTPUT), split = "\\ "), "[[",19))

    #if(grepl(INVAR_score.path, pattern = "downsampled") == TRUE){
    #    print("Downsampled molecules data (downsampled to X GE)")
    # the downsampled # mutations columns was used to state target GE
    #  colnames(INVAR_scores)[which(colnames(INVAR_scores) == "n_mutations")] <- "target_GE"
    #  INVAR_scores$iteration <- unlist(lapply(strsplit(as.character(INVAR_scores$OUTPUT), split = "\\ "), "[[",19))


  } else {
    print("wrong column length")
    return(NULL)
  }

  INVAR_scores$overall_filter <- paste0(INVAR_scores$GLRT_setting, "_", INVAR_scores$POLISH.PASS, "_", INVAR_scores$BOTH_STRANDS.PASS, "_", INVAR_scores$LOCUS_NOISE.PASS)

  print("annotating with SLX_barcode")
  INVAR_scores$SLX_barcode <- gsub(INVAR_scores$FILE_NAME, pattern = SLX_barcode_format, replacement = "\\1")
  INVAR_scores$SLX_barcode <- unlist(lapply(strsplit(INVAR_scores$SLX_barcode, split = "\\."), "[[", 1))
  INVAR_scores$data <- gsub(INVAR_scores$FILE_NAME, pattern = ".*\\.(.*ptspec)\\..*", replacement = "\\1")

  # get initial rows
  initial_rows <- nrow(INVAR_scores)

  ## add SLX_barcode if it's not present
  if(length(which(colnames(SLX_layout) %in% "SLX_barcode")) == 1){
    print("SLX_barcode is present, good")
  } else{
    print("SLX_barcode column is not present, adding it now")
    SLX_layout$barcode <- gsub(SLX_layout$barcode, pattern = "-", replacement = "_")
    SLX_layout$SLX_barcode <- paste0(SLX_layout$SLX_ID, "_", SLX_layout$barcode)
  }

  ## merge
  print("merging")
  INVAR_scores <- left_join(INVAR_scores, SLX_layout[,c("Study", "SLX_barcode", "case_or_control", "sample_type", "data_type", "Patient", "Timepoint")], by = "SLX_barcode")

  INVAR_scores <- dplyr::select(INVAR_scores, - FILE_NAME, - OUTPUT)
  INVAR_scores$adjusted_INVAR <- INVAR_scores$INVAR_SCORE
  INVAR_scores$adjusted_IMAF <- INVAR_scores$IMAF


  if(adjust == TRUE){
    INVAR_scores[INVAR_scores$mutant_present == FALSE, "adjusted_INVAR"] <- 0
    INVAR_scores[INVAR_scores$mutant_present == FALSE, "adjusted_IMAF"] <- 0

    print("setting ctDNA level to 1/# molecules where p_mle < 1/# molecules")
    INVAR_scores[INVAR_scores$adjusted_IMAF < 1/INVAR_scores$DP & INVAR_scores$adjusted_IMAF != 0,'adjusted_IMAF'] <-
      1/INVAR_scores[INVAR_scores$adjusted_IMAF < 1/INVAR_scores$DP & INVAR_scores$adjusted_IMAF != 0,'DP']
    #print("ctDNA level = average ")
    #INVAR_scores$adjusted_IMAF <- INVAR_scores$mut_sum/ INVAR_scores$DP
  }

  ptspec <- filter(INVAR_scores, data == "ptspec")
  nonptspec <- filter(INVAR_scores, data == "nonptspec")


  INVAR_scores <- rbind(ptspec, nonptspec)

  INVAR_scores.conditions <- unique(INVAR_scores$overall_filter)

  INVAR_score <- data.frame()
  threshold_df <- data.frame()

  INVAR_score.raw <<- INVAR_scores

  for (k in 1:length(INVAR_scores.conditions)) {

    curr <- filter(INVAR_scores, overall_filter == INVAR_scores.conditions[k] & data == "ptspec" | overall_filter == INVAR_scores.conditions[k] & data == "nonptspec" & CONTAMINATION_RISK.PASS == TRUE)

    if(adjust == TRUE){
      print("using adjusted invar score")
      threshold_curr <- filter(curr, data == "nonptspec")$adjusted_INVAR[order(-filter(curr, data == "nonptspec")$adjusted_INVAR)][ceiling((1-INVAR_score_specificity)*nrow(filter(curr, data == "nonptspec")))]
    } else{
      print("not using adjusted invar score")
      threshold_curr <- filter(curr, data == "nonptspec")$INVAR_SCORE[order(-filter(curr, data == "nonptspec")$INVAR_SCORE)][ceiling((1-INVAR_score_specificity)*nrow(filter(curr, data == "nonptspec")))]
    }

    print(paste("Threshold for", INVAR_scores.conditions[k], "is" , threshold_curr))

    threshold_df <- rbind(data.frame(INVAR_scores.conditions[k], threshold_curr), threshold_df)

    curr$detection <- curr$adjusted_INVAR > threshold_curr
    INVAR_score <- rbind(INVAR_score, curr)
  }
  threshold_df  <- threshold_df

  INVAR_score$data_code <- 0
  INVAR_score$data_code[INVAR_score$data == "ptspec"] <- 1

  # additional variables
  INVAR_score$sample_name <- paste0(INVAR_score$patient, " (", INVAR_score$pt_mut_belongs_to, ")")
  INVAR_score$additional_filter <- paste0(INVAR_score$POLISH.PASS, "_", INVAR_score$BOTH_STRANDS.PASS, "_", INVAR_score$LOCUS_NOISE.PASS, "_", INVAR_score$data)

  # check rows
  final_rows <- nrow(INVAR_scores)
  print(paste(initial_rows, final_rows))

  if (!is.na(filter.pattern)){
  INVAR_score <- INVAR_score[grepl(INVAR_score$overall_filter, pattern = filter.pattern),]
  }

   return(INVAR_score)
}

reformat.GLRT.comparison <- function(INVAR_score){
  print("reformating GLRT settings to compare between using size and not using size")
  INVAR_score.size <- filter(INVAR_score, GLRT_setting == "INVAR_SCORE.using_size")
  INVAR_score.nosize <- filter(INVAR_score, GLRT_setting == "INVAR_SCORE.no_size")

  INVAR_score.GLRT.comparison <- left_join(INVAR_score.size, INVAR_score.nosize, by = c("patient", "pt_mut_belongs_to", "POLISH.PASS", "BOTH_STRANDS.PASS", "LOCUS_NOISE.PASS", "smooth", "called_against", "CONTAMINATION_RISK.PASS", "outlier_suppression", "null_likelihood", "alt_likelihood", "SLX_barcode", "data", "Study", "sample_name", "case_or_control", "sample_type", "data_type", "Patient", "Timepoint","data_code", "additional_filter"))

  INVAR_score.GLRT.comparison$detection_change <- paste0(INVAR_score.GLRT.comparison$detection.x, "_", INVAR_score.GLRT.comparison$detection.y)

  detection_change_names <<- c("FALSE_FALSE" = "Never detected",
                               "FALSE_TRUE" = "Detected without size",
                               "TRUE_FALSE" = "Detected with size",
                               "TRUE_TRUE" = "Always detected")

  return(INVAR_score.GLRT.comparison)
}

# Downsampling of GE functions ------------------------
downsample_loci_fn <- function(data.pre_downsample, x){
  # dilution series of MR1022
  # BRAF should not be locus noise fail. this is the only exception because it is ubiquitous
  # x = number of mutations
  library(dplyr)
  source("~/tapas_pipeline/R/TAPAS_functions.R")

  print("setting BRAF to locus noise pass")
  data.pre_downsample[data.pre_downsample$uniq_pos == "chr7:140453136", 'LOCUS_NOISE.PASS'] <- T
  data <- filter(data.pre_downsample, LOCUS_NOISE.PASS == T)

  #("ensuring all loci are in the mutations csv")
  tumour_mutations <- read.csv(file = "/scratchb/nrlab/wan01/JXPT/JXP0268/bed/MELR.patient_specific.baseline.csv")
  print(paste("tumour mutations rows", nrow(tumour_mutations)) )

  tumour_mutations <- filter(tumour_mutations, Patient == "MR1022")
  #print(paste("tumour mutations rows", nrow(tumour_mutations)))
  tumour_mutations <- tumour_mutations[tumour_mutations$uniq_pos %in% data$uniq_pos,]
  #print(paste("tumour mutations rows", nrow(tumour_mutations)))
  tumour_AF_key <- tumour_mutations[order(-tumour_mutations$tumour_AF),]

  print(paste("taking the top", x , "loci for this mixture sample"))
  curr_loci <- tumour_AF_key[1:x, 'uniq_pos']
  data <- data[data$uniq_pos %in% curr_loci ,]

  return(data)
}
# downsample families to a fixed number of final GE
downsample_families <- function(raw, final_GE = 10){

  # represent mutant reads only once
  represent_mut_reads_once_only <- function(raw.pass){
    print("mutant reads in the filtered data frame will only be observed once for this downsampling analysis")
    raw.pass.mutant <- filter(raw.pass, mutant == 1) %>% distinct()
    raw.pass.wt <- filter(raw.pass, mutant == 0)

    raw.pass <- rbind(raw.pass.wt, raw.pass.mutant)
    return(raw.pass)
  }

  raw <- represent_mut_reads_once_only(raw)

  # N.B. effectively 10 available GE per droplet of blood
  median_DP <- median(data.frame(table(raw$uniq_pos))[,2])

  downsample_factor <- median_DP/final_GE
  total_dp_downsampled <- round(nrow(raw)/downsample_factor)


  downsampled <- raw[sample(nrow(raw), total_dp_downsampled, replace = T),]
  table(downsampled$mutant)

  return(downsampled)
}

# downsample families to a % of the starting GE
downsample_families.proportion <- function(raw, proportion){

  total_dp_downsampled <- round(nrow(raw) * proportion)

  downsampled <- raw[sample(nrow(raw), total_dp_downsampled, replace = T),]

  #median(data.frame(table(downsampled$uniq_pos))[,2])
  return(downsampled)
}
