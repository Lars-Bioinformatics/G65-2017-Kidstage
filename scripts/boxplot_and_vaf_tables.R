library(tidyverse)

save.image("debug.RData")

# data = read.table("/work/G65-2017-Kidstage/Connor/mutect2_somatic_variants/bam_readcount_output_q20_b20/G65-T10-BA34_nimblegen-medexome_HYLKFDSXX/combined_G65-T10-BA34_nimblegen-medexome_HYLKFDSXX_bam-readcount_VAF_q20_b20.txt", header=T)
# print(snakemake)
print(snakemake@input$vaf)
data = do.call(rbind, lapply(snakemake@input$vaf, read.table, header=T))
print(table(data$tumor))

print(snakemake@output$table_long)
print(snakemake@output)

# Summarised table of VAF (long format)
dt_long = data %>% 
    group_by(tumor, tissue, chr, pos, ref, alt) %>%
    summarise(mean_RD=mean(RD), mean_AD=mean(AD), mean_VAF=mean(VAF), .groups="drop")
write.table(dt_long, snakemake@output$table_long, row.names=F, quote=F)

# In wide format
dt_wide = pivot_wider(dt_long, names_from=tissue, values_from=c(mean_RD,mean_AD,mean_VAF))
write.table(dt_wide, snakemake@output$table_wide, row.names=F, quote=F)

# boxplot
p = data %>% ggplot() +
    geom_boxplot(aes(tumor,VAF,color=tissue)) +
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    labs(title = "Tumor Signal in Plasma",
         x="Sample",
         y="VAF")
p
ggsave(snakemake@output$boxplot,width = 32, height = 9)