# Load packages, and install if not
packages <- c("stringr", "rebus", "yaml")

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE, repos = c("https://mirrors.dotsrc.org/cran/", "http://cran.us.r-project.org"))
  require(p, character.only = TRUE)
}

for (i in 1:length(packages)) {
  usePackage(packages[i])
}

# Read CLI args
yaml <- read_yaml("config.yaml")

# List files
files <- list.files(yaml$inputdir, "*.bam")
files <- str_subset(files, ".bam$")
files <- sub(yaml$samplesubset$bamsuffix, "", files)

# Divide into tumor and refsample samples
refsample <- str_subset(files, yaml$samplesubset$refsample)
tumor <- str_subset(files, yaml$samplesubset$tumor)
plasma <- str_subset(tumor, yaml$samplesubset$plasma)
biopsy <- str_subset(tumor, yaml$samplesubset$biopsy)

# Create regex pattern
pattern <- START %R%
  one_or_more(char_class(ASCII_ALNUM)) %R%
  "-" %R%
  one_or_more(DGT) %R%
  "-"

# Matched samples
func_matched <- function(vec_samples, tsv_out) {
  tmp_matrix <- matrix(nrow = length(vec_samples),
                       ncol = 4)
  tmp_matrix[1:length(vec_samples),1] <- vec_samples
  for (i in 1:length(vec_samples)) {
    tmp_matrix[i,2] <- refsample[str_detect(refsample, str_extract(tmp_matrix[i,1], pattern))][1]
    tmp_matrix[i,3] <- gsub('.{1}$', '', str_extract(tmp_matrix[i,1], pattern))
    tmp_matrix[i,4] <- str_replace(tmp_matrix[i,1], pattern, "-")
  }
  write.table(tmp_matrix, 
              file = tsv_out, 
              sep = "\t", 
              quote = F,
              row.names = F,
              col.names = F)
}
func_matched(tumor, "config/matched-tumor-normal.tsv")
func_matched(plasma, "config/matched-plasma-normal.tsv")
func_matched(biopsy, "config/matched-biopsy-normal.tsv")

# Ensemble data frame
func_ensemble <- function(vec_samples, tsv_out) {
  VARIANTCALLEROUTPUTS <- list()
  VARIANTCALLERS <- list()
  for (i in 1:length(yaml$variantcaller)) {
    tmp <- unlist(yaml$variantcaller[i])
    if (as.logical(tmp[1])){
      if (tmp[2] == "snvs") {
        VARIANTCALLEROUTPUTS[length(VARIANTCALLEROUTPUTS)+1] <- paste0(str_replace(names(tmp)[1], ".include", ""), "_snvs.vcf")
      } else if (tmp[2] == "indels") {
        VARIANTCALLEROUTPUTS[length(VARIANTCALLEROUTPUTS)+1] <- paste0(str_replace(names(tmp)[1], ".include", ""), "_indels.vcf")
      } else if (tmp[2] == "both") {
        VARIANTCALLEROUTPUTS[length(VARIANTCALLEROUTPUTS)+1] <- paste0(str_replace(names(tmp)[1], ".include", ""), "_snvs.vcf")
        VARIANTCALLEROUTPUTS[length(VARIANTCALLEROUTPUTS)+1] <- paste0(str_replace(names(tmp)[1], ".include", ""), "_indels.vcf")
      } else if (tmp[2] == "bothcombined") {
        VARIANTCALLEROUTPUTS[length(VARIANTCALLEROUTPUTS)+1] <- paste0(str_replace(names(tmp)[1], ".include", ""), "_bothcombined.vcf")
      }
    }
    if (as.logical(tmp[1])){
      VARIANTCALLERS[length(VARIANTCALLERS)+1] <- str_replace(names(tmp)[1], ".include", "")
    }
  }
  VARIANTCALLEROUTPUTS <- unlist(VARIANTCALLEROUTPUTS)
  VARIANTCALLERS <- unlist(VARIANTCALLERS)
  tmp_matrix <- matrix(nrow = length(vec_samples),
                       ncol = 4)
  tmp_matrix[1:length(vec_samples),1] <- vec_samples
  for (i in 1:length(vec_samples)) {
    tmp_matrix[i,2] <- refsample[str_detect(refsample, str_extract(tmp_matrix[i,1], pattern))][1]
    tmp_matrix[i,3] <- gsub('.{1}$', '', str_extract(tmp_matrix[i,1], pattern))
    tmp_matrix[i,4] <- str_replace(tmp_matrix[i,1], pattern, "-")
  }
  tmp_matrix <- do.call(rbind, replicate(length(VARIANTCALLEROUTPUTS), tmp_matrix, simplify = F))
  tmp_matrix <- cbind(tmp_matrix, 
                      rep(unlist(VARIANTCALLEROUTPUTS), each = length(vec_samples)),
                      rep(unlist(VARIANTCALLEROUTPUTS), each = length(vec_samples)))
  for (i in 1:nrow(tmp_matrix)) {
    tmp_matrix[i,6] <- paste0(tmp_matrix[i,5], tmp_matrix[i,3])
  }
  write.table(tmp_matrix, 
              file = tsv_out, 
              sep = "\t", 
              quote = F,
              row.names = F,
              col.names = F)
}
func_ensemble(biopsy, "config/ensemble-biopsy-normal.tsv")

