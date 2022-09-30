#===================================================================================
#
# FILE: config2.R
# USAGE: Please update variables in this R script for filtering of INVAR data points.
#        Default values are given in square brackets
#
#===================================================================================

# User to define these variables below ========================
# Filtering INVAR loci
max_DP <- 2000 #[1500] Omit data points with uncharacteristially high unique DP given the input mass used
min_ref_DP <- 10 #[5] Min DP is set by mpileups as 5, here we require at least 5 ref reads at a locus. Change to 0 for sWGS data
individual_MQSB_threshold <-  0.01 #[0.01] Excludes data points due to poor MQ and SB, but locus is retained
n_alt_alleles_threshold <-  3 #[3] Blacklist loci with >= N separate alternate alleles.
minor_alt_allele_threshold <- 2 #[2] Blacklist multiallelic loci with a mutant read count of >= N in the minor mutant allele.
FR_ratio.min_reads    <- 5 #[10] Require a minimum number of mutant reads at a locus in order to call it as strand-biased.
FR_ratio.threshold    <- 0.2 #[0.2] The ratio between F and R mutant reads required to indicate that there is bias between F and R
proportion_of_controls <- 0.1 #[0.1] Blacklist loci that have signal in >30% of the nonptspec samples
max_background_mean_AF <-  0.01 #[0.01] filter loci with a background AF in controls greater than this value, change to 0.1 for sWGS samples
cosmic_threshold <- 0 # [0] Loci with >0 entries in COSMIC are considered as COSMIC mutations
outlier_suppression_threshold_settings <- 0.05 # this can be a vector if you wish to vary its stringency.
only_weigh_mutants <- T # [TRUE] Only weigh ctDNA signal based on mutant fragments
is.blood_spot <- F # [FALSE] Only change to TRUE if you are running blood spot or sWGS data through the pipeline. This omits outlier-suppression on samples with deduplicated depth of <5x because high AF loci cannot be reliably identified with low depth

# handling different SLX barcode formats
barcode_format <- "SXTHS[A-Z]0[0-9]" # "D.*_D.*\\." the user needs to define the barcode structure here in a regular expression
SLX_format <- "SLX-[0-9]+"
SLX_barcode_format <- paste0(".*(", SLX_format, "_" , barcode_format, ").*")

# Once the Rdata file has been generated it can be loaded directly. Set this parameter below
load_existing_Rdata <- F

# File paths and suffixes (no need to edit unless altering file structure in repo)========================
input_dir <- "input/"
gz_dir <- "output_gz/"
output_dir <- "output_R/"
TAPAS.file_suffix <- ".combined.final.ann.tsv.gz"
TAPAS.tsv_colnames <- c('CHROM', 'POS', 'REF', 'ALT', 'DP', 'DP4','REF_F', 'ALT_F','REF_R', 'ALT_R', 'MQSB', 'SLX', 'BARCODE', 'COSMIC_MUTATIONS', 'COSMIC_SNP', '1KG_AF', 'TRINUCLEOTIDE', 'FILE_NAME')
