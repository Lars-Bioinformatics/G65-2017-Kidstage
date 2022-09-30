# Load packages, and install if not
packages <- c("data.table", "rebus", "dplyr", "plyr", "ggplot2", "stringr", "stringi", "devtools", "grid", "gridExtra")

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE, repos = c("https://mirrors.dotsrc.org/cran/", "http://cran.us.r-project.org"), INSTALL_opts = c('--no-lock'))
  require(p, character.only = TRUE)
}

for (i in 1:length(packages)) {
  usePackage(packages[i])
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", INSTALL_opts = c('--no-lock'))
if (!requireNamespace("VariantAnnotation", quietly = TRUE))
  BiocManager::install("VariantAnnotation", INSTALL_opts = c('--no-lock'))
library("VariantAnnotation")

# Import files
args = commandArgs(trailingOnly=TRUE)

pattern_extension <- DOT %R%
  one_or_more(char_class(ASCII_ALNUM)) %R%
  END
samples <- str_remove(basename(args[2:length(args)]), pattern_extension)

pattern_patients <- START %R%
  one_or_more(char_class(ASCII_ALNUM)) %R%
  "-" %R%
  one_or_more(DGT) %R%
  "-"
patients <- unique(str_extract(samples, pattern_patients))

tsv_list <- list()
filesin <- list()
case <- list()
case_samples <- list()
case_zscore <- list()
control <- list()
control_sd <- list()
control_mean <- list()
sites_detected <- list()
for (i in 1:(length(args)-1)) {
  tsv_list[[i]] <- list()
  sites_detected[[i]] <- list()
  case_zscore[[i]] <- list()
  filesin[[i]] <- list.files(path = paste0(args[i+1], "/"), pattern = DOT %R% "vcf" %R% END)
  filesin[[i]] <- filesin[[i]][!str_detect(filesin[[i]], "qual")]
  case[[i]] <- str_detect(filesin[[i]], patients[i])
  case_samples[[i]] <- filesin[[i]][case[[i]]]
  control[[i]] <- !str_detect(filesin[[i]], patients[i])
  for (j in 1:length(filesin[[i]])){
    tsv_list[[i]][[j]] <- VariantAnnotation::readVcf(file = paste0(paste0(args[i+1], "/"), filesin[[i]][j]), genome = "hg38", param = ScanVcfParam())
    sites_detected[[i]][[j]] <- mean(unlist(assays(tsv_list[[i]][[j]])$AF[,2]))
    sites_detected[[i]] <- unlist(sites_detected[[i]])
  }
  control_sd[[i]] <- sd(sites_detected[[i]][control[[i]]])
  control_mean[[i]] <- mean(sites_detected[[i]][control[[i]]])
  for (j in 1:length(filesin[[i]][case[[i]]])) {
    case_zscore[[i]][[j]] <- (sites_detected[[i]][case[[i]]][j] - control_mean[[i]]) / control_sd[[i]]
  }
}

pattern_samples <- START %R%
  one_or_more(char_class(ASCII_ALNUM)) %R%
  "-" %R%
  one_or_more(DGT) %R%
  "-" %R%
  "plasma" %R%
  one_or_more(ascii_digit())

pattern_method <- capture(one_or_more(char_class(ASCII_ALNUM) %R% "_-"))

pattern_method_sample <- "/" %R%
  pattern_method %R%
  "/" %R%
  one_or_more(char_class(ASCII_ALNUM %R% "_-")) %R%
  optional("/") %R%
  END

samples <- str_match(unlist(case_samples), pattern_samples)

snv_method <- stri_match_first_regex(args[2], pattern_method_sample)[1,2]

pattern_patients_date <- START %R% 
  capture(one_or_more(char_class(ASCII_ALNUM)) %R%
            "-" %R%
            one_or_more(DGT)) %R%
  "-" %R%
  one_or_more(char_class(ASCII_ALNUM)) %R%
  capture(repeated(DGT, 2, 2)) %R%
  capture(repeated(DGT, 2, 2)) %R%
  capture(repeated(DGT, 2, 2))

patients_date <- stri_match_first_regex(samples, pattern_patients_date)

out_df <- data.frame("patients" = patients_date[,2],
                     "samples" = paste0(patients_date[,2], "  20", patients_date[,3], "-", patients_date[,4], "-", patients_date[,5]),
                     "Zscores" = unlist(case_zscore))

write.csv(out_df, paste0(args[1], "snv_", snv_method, ".csv"))

min_zscore <- 1.645 # probability of 0.95

p <- ggplot(data = out_df, aes(x = samples,
                               y = Zscores,
                               fill = patients)) +
  xlab("Date") +
  ylab("Z-score") +
  geom_bar(position="dodge", 
           stat="identity",
           color = "black") +
  geom_hline(yintercept = min_zscore,
             lty = 2) +
  scale_fill_discrete("Patients") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(paste0(args[1], "snv_", snv_method, ".pdf"),
       p,
       device = pdf(),
       width = 5,
       height = 4)

dev.off()


tsv_list <- list()
filesin <- list()
case <- list()
case_samples <- list()
case_zscore <- list()
control <- list()
control_sd <- list()
control_mean <- list()
sites_detected <- list()
for (i in 1:(length(args)-1)) {
  tsv_list[[i]] <- list()
  sites_detected[[i]] <- list()
  case_zscore[[i]] <- list()
  filesin[[i]] <- list.files(path = paste0(args[i+1], "/"), pattern = DOT %R% "vcf" %R% END)
  filesin[[i]] <- filesin[[i]][str_detect(filesin[[i]], "qual")]
  case[[i]] <- str_detect(filesin[[i]], patients[i])
  case_samples[[i]] <- filesin[[i]][case[[i]]]
  control[[i]] <- !str_detect(filesin[[i]], patients[i])
  for (j in 1:length(filesin[[i]])){
    tsv_list[[i]][[j]] <- VariantAnnotation::readVcf(file = paste0(paste0(args[i+1], "/"), filesin[[i]][j]), genome = "hg38", param = ScanVcfParam())
    sites_detected[[i]][[j]] <- mean(unlist(assays(tsv_list[[i]][[j]])$AF[,2]))
    sites_detected[[i]] <- unlist(sites_detected[[i]])
  }
  control_sd[[i]] <- sd(sites_detected[[i]][control[[i]]])
  control_mean[[i]] <- mean(sites_detected[[i]][control[[i]]])
  for (j in 1:length(filesin[[i]][case[[i]]])) {
    case_zscore[[i]][[j]] <- (sites_detected[[i]][case[[i]]][j] - control_mean[[i]]) / control_sd[[i]]
  }
}

pattern_samples <- START %R%
  one_or_more(char_class(ASCII_ALNUM)) %R%
  "-" %R%
  one_or_more(DGT) %R%
  "-" %R%
  "plasma" %R%
  one_or_more(ascii_digit())

pattern_method <- capture(one_or_more(char_class(ASCII_ALNUM) %R% "_-"))

pattern_method_sample <- "/" %R%
  pattern_method %R%
  "/" %R%
  one_or_more(char_class(ASCII_ALNUM %R% "_-")) %R%
  optional("/") %R%
  END

samples <- str_match(unlist(case_samples), pattern_samples)

snv_method <- stri_match_first_regex(args[2], pattern_method_sample)[1,2]

pattern_patients_date <- START %R% 
  capture(one_or_more(char_class(ASCII_ALNUM)) %R%
            "-" %R%
            one_or_more(DGT)) %R%
  "-" %R%
  one_or_more(char_class(ASCII_ALNUM)) %R%
  capture(repeated(DGT, 2, 2)) %R%
  capture(repeated(DGT, 2, 2)) %R%
  capture(repeated(DGT, 2, 2))

patients_date <- stri_match_first_regex(samples, pattern_patients_date)

out_df <- data.frame("patients" = patients_date[,2],
                     "samples" = paste0(patients_date[,2], "  20", patients_date[,3], "-", patients_date[,4], "-", patients_date[,5]),
                     "Zscores" = unlist(case_zscore))

write.csv(out_df, paste0(args[1], "qual_snv_", snv_method, ".csv"))

min_zscore <- 1.645 # probability of 0.95

p <- ggplot(data = out_df, aes(x = samples,
                               y = Zscores,
                               fill = patients)) +
  xlab("Date") +
  ylab("Z-score") +
  geom_bar(position="dodge", 
           stat="identity",
           color = "black") +
  geom_hline(yintercept = min_zscore,
             lty = 2) +
  scale_fill_discrete("Patients") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(paste0(args[1], "qual_snv_", snv_method, ".pdf"),
       p,
       device = pdf(),
       width = 5,
       height = 4)

dev.off()

