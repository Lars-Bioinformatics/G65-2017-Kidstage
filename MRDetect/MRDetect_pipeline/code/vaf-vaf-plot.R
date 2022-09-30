# Load packages, and install if not
packages <- c("data.table", "rebus", "dplyr", "plyr", "ggplot2", "stringr")

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE, repos = c("https://mirrors.dotsrc.org/cran/", "http://cran.us.r-project.org"), INSTALL_opts = c('--no-lock'))
  require(p, character.only = TRUE)
}

for (i in 1:length(packages)) {
  usePackage(packages[i])
}

# Import files
args = commandArgs(trailingOnly=TRUE)

filepath <- ifelse(str_sub(args[1], -1) == "/",
                   args[1],
                   paste0(args[1], "/"))

tsv_list = list()
for (i in 1:(length(args)-1)) {
  tsv_list[[i]] <- fread(args[i+1])
  tsv_list[[i]] <- tsv_list[[i]][tsv_list[[i]]$VCF_SAMPLE_ID == "TUMOR",]
}

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

vafvafplot <- function(df_1, df_2, samplename_1, samplename_2) {
  # Collapse CHROM and POS of df_1
  tmp_1 <- as.matrix(rep(NA, nrow(df_1)))
  for (i in 1:nrow(df_1)) {
    tmp_1[i, 1] <- paste0(df_1[i, 1:2], collapse = "-")
  }
  tmp_1 <- data.frame("CHROM.POS" = tmp_1,
                      "VAF" = df_1$AF)
  # Collapse CHROM and POS of df_2
  tmp_2 <- as.matrix(rep(NA, nrow(df_2)))
  for (i in 1:nrow(df_2)) {
    tmp_2[i, 1] <- paste0(df_2[i, 1:2], collapse = "-")
  }
  tmp_2 <- data.frame("CHROM.POS" = tmp_2,
                      "VAF" = df_2$AF)
  # Create combined df
  tmp_df <- full_join(tmp_1, 
                      tmp_2, 
                      by = "CHROM.POS",
                      suffix = c("1", "2"))
  tmp_df[,2] <- as.numeric(tmp_df[,2])
  tmp_df[,3] <- as.numeric(tmp_df[,3])
  ggplot(tmp_df, aes(x = VAF1,
                     y = VAF2)) +
    geom_point(color = alpha("deepskyblue2", 0.5)) +
    geom_abline(slope = 1,
                intercept = 0,
                lwd = 0.25,
                lty = 3) +
    labs(x = bquote("VAF"[.(samplename_1)] ~ " ( n =" ~ .(length(tmp_df[,2][tmp_df[,2] != 0])) ~ ")"),
         y = bquote("VAF"[.(samplename_2)] ~ " ( n =" ~ .(length(tmp_df[,3][tmp_df[,3] != 0])) ~ ")")) +
    theme_bw()
}

for (i in 1:length(patients)) {
  outdir <- paste0(filepath, gsub('.{1}$', '', patients[i]), "/")
  if (!dir.exists(paste0(outdir))) {
    dir.create(outdir)
  }
  samples_indices <- str_which(samples, pattern = patients[i])
  samples_list <- list()
  for (j in 1:length(samples_indices)) {
    samples_list[[j]] <- tsv_list[[samples_indices[j]]]
  }
  samples_vector <- str_subset(samples, pattern = patients[i])
  j = 1; k = 2; while(j < length(samples_list)) {
    p <- vafvafplot(samples_list[[j]], samples_list[[k]], samples_vector[j], samples_vector[k])
    ggsave(filename = paste0(samples_vector[j], "_vs_", samples_vector[k], ".pdf"),
           plot = p,
           device = pdf(),
           path = outdir,
           width = 12,
           height = 12,
           units = "cm")
    if (k == length(samples_list)) {
      j = j + 1
      k = j + 1
    } else {
      k = k + 1
    }
    dev.off()
  }
}
