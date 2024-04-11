#---- Load libraries------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggrepel)
#---- Parameters ---------------------------------------------------------------
cutoff_FC <- 2
cutoff_FDR <- 0.05
#---- Result from DESeq2 -------------------------------------------------------
res <- read.csv("differential_analysis.csv")
#---- Label top 10 genes -------------------------------------------------------
res <- res %>%
  mutate(label = if_else(row_number() <= 10, as.character(gene_name), ""))
#---- Volcano plot -------------------------------------------------------------
ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = significant, label = label)) +
  xlab("log2(FC)") + 
  ylab("-log10(FDR)") +
  scale_color_manual(values = c("blue3", "grey", "red3")) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = log2(cutoff_FC), linetype=3, colour = "grey60") +
  geom_vline(xintercept = -log2(cutoff_FC), linetype=3, colour = "grey60") +
  geom_hline(yintercept = -log10(cutoff_FDR), linetype = 3, colour = "grey60") +
  geom_point(show.legend = FALSE) +
  geom_label_repel(force = 1, size = 3, box.padding = 0.5, 
                   point.padding = 0.3, show.legend=FALSE) +
  theme_bw()
ggsave('volcano_plot.png', width = 5, height = 5)

