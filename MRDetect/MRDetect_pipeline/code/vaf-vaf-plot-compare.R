# Load packages, and install if not
packages <- c("data.table", "rebus", "dplyr", "plyr", "ggplot2", "stringr", "devtools", "grid", "gridExtra")

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

list_tumor = list()
list_normal = list()
for (i in 1:(length(args)-1)) {
  tmp <- VariantAnnotation::readVcf(file = args[i+1], genome = "hg38", param = ScanVcfParam())
  list_tumor[[i]] <- data.frame("ID" = names(assays(tmp)$VAF[,2]), 
                                "VAF" = assays(tmp)$VAF[,2]) 
  list_normal[[i]] <- data.frame("ID" = names(assays(tmp)$VAF[,1]), 
                                 "VAF" = assays(tmp)$VAF[,1]) 
}

pattern_extension <- DOT %R%
  one_or_more(char_class(ASCII_ALNUM)) %R%
  END
samples <- str_remove(basename(args[2:length(args)]), pattern_extension)

vafvafplot <- function(df_1_1, df_1_2, df_2_1, df_2_2, samplename_1, samplename_2) {
  # Create combined df
  tmp_df <- full_join(df_1_1, 
                      df_2_1, 
                      by = "ID",
                      suffix = c("1", "2"))
  tmp_df[,2] <- as.numeric(tmp_df[,2])
  tmp_df[,3] <- as.numeric(tmp_df[,3])
  for (i in 1:nrow(tmp_df)) {
    if (is.na(tmp_df[i,2]))
      tmp_df[i,2] <- 0
    if (is.na(tmp_df[i,3]))
      tmp_df[i,3] <- 0
  }
  booldf1 <- rep(FALSE, nrow(df_1_1))
  for (i in 1:nrow(df_1_1)) {
    if (dim(df_1_1[df_2_1$ID == df_1_1$ID[i],])[1] == 0)
      booldf1[i] <- TRUE
  }
  booldf2 <- rep(FALSE, nrow(df_2_1))
  for (i in 1:nrow(df_2_1)) {
    if (dim(df_2_1[df_1_1$ID == df_2_1$ID[i],])[1] == 0)
      booldf2[i] <- TRUE
  }
  p1 <- ggplot(tmp_df, aes(x = VAF1, y = VAF2)) +
    geom_point(color = alpha("deepskyblue2", 0.5)) +
    geom_abline(slope = 1,
                intercept = 0,
                lwd = 0.25,
                lty = 3) +
    xlim(-0.01, max(c(df_1_1$VAF, df_2_1$VAF, df_1_2$VAF, df_2_2$VAF))+0.01) + 
    ylim(-0.01, max(c(df_1_1$VAF, df_2_1$VAF, df_1_2$VAF, df_2_2$VAF))+0.01) +
    labs(x = bquote("VAF"[.(samplename_1)] ~ " ( n =" ~ .(length(tmp_df[,2][tmp_df[,2] != 0])) ~ ")"),
         y = bquote("VAF"[.(samplename_2)] ~ " ( n =" ~ .(length(tmp_df[,3][tmp_df[,3] != 0])) ~ ")")) +
    theme_bw()
  
  p2_1 <- ggplot(df_1_1, aes(x = VAF)) +
    geom_histogram(color = "black",
                   fill = alpha("deepskyblue2"),
                   binwidth = 0.01) +
    geom_histogram(data = df_1_1[booldf1,],
                   color = "black",
                   fill = alpha("orchid2"),
                   binwidth = 0.01) +
    xlim(-0.01, max(c(df_1_1$VAF, df_2_1$VAF, df_1_2$VAF, df_2_2$VAF))+0.01) + 
    labs(x = bquote("VAF"[.(samplename_1) ~ "-Tumor"]),
         y = bquote("Count"[" "])) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  
  p2_2 <- ggplot(df_1_2, aes(x = VAF)) +
    geom_histogram(color = "black",
                   fill = alpha("deepskyblue2"),
                   binwidth = 0.01) +
    geom_histogram(data = df_1_2[booldf1,],
                   color = "black",
                   fill = alpha("orchid2"),
                   binwidth = 0.01) +
    xlim(-0.01, max(c(df_1_1$VAF, df_2_1$VAF, df_1_2$VAF, df_2_2$VAF))+0.01) +
    labs(x = bquote("VAF"[.(samplename_1) ~ "-Normal"]),
         y = bquote("Count"[" "])) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  
  p3_1 <- ggplot(df_2_1, aes(x = VAF)) +
    geom_histogram(color = "black",
                   fill = alpha("deepskyblue2"),
                   binwidth = 0.01) +
    geom_histogram(data = df_2_1[booldf2,],
                   color = "black",
                   fill = alpha("orchid2"),
                   binwidth = 0.01) +
    xlim(-0.01, max(c(df_1_1$VAF, df_2_1$VAF, df_1_2$VAF, df_2_2$VAF))+0.01) +
    labs(x = bquote("VAF"[.(samplename_2) ~ "-Tumor"]),
         y = bquote("Count"[" "])) +
    theme_bw() +
    theme(panel.grid.major.y  = element_blank(), 
          panel.grid.minor.y = element_blank()) +
    coord_flip()
  
  p3_2 <- ggplot(df_2_2, aes(x = VAF)) +
    geom_histogram(color = "black",
                   fill = alpha("deepskyblue2"),
                   binwidth = 0.01) +
    geom_histogram(data = df_2_2[booldf2,],
                   color = "black",
                   fill = alpha("orchid2"),
                   binwidth = 0.01) +
    xlim(-0.01, max(c(df_1_1$VAF, df_2_1$VAF, df_1_2$VAF, df_2_2$VAF))+0.01) +
    labs(x = bquote("VAF"[.(samplename_2) ~ "-Normal"]),
         y = bquote("Count"[" "])) +
    theme_bw() +
    theme(panel.grid.major.y  = element_blank(), 
          panel.grid.minor.y = element_blank()) +
    coord_flip()
  
  layout_grid <- rbind(c(2,2,2,1,1),
                       c(3,3,3,1,1),
                       c(4,4,4,5,6),
                       c(4,4,4,5,6),
                       c(4,4,4,5,6))
  grid.arrange(grid.rect(gp=gpar(col="white")), p2_2, p2_1, p1, p3_1, p3_2, layout_matrix = layout_grid)
}

j = 1; k = 2; while(j < length(list_tumor)) {
  pdf(file = paste0(args[1], samples[j], "_vs_", samples[k], ".pdf"),
      width = 10,
      height = 10)
  vafvafplot(list_tumor[[j]], list_normal[[j]], list_tumor[[k]], list_normal[[k]], samples[j], samples[k])
  dev.off()
  if (k == length(list_tumor)) {
    j = j + 1
    k = j + 1
  } else {
    k = k + 1
  }
}
