#---- Load libraries -----------------------------------------------------------
library(DESeq2)
library(edgeR)
library(dplyr)
library(ggplot2)
#---- Input and output directories ---------------------------------------------
# Define directory paths
dir_data <- './data/'
dir_output <- './output/'
# Check if the output directory exists, if not create it
if (!dir.exists(dir_output)) {
  dir.create(dir_output)
}
#---- Raw data -----------------------------------------------------------------
sampleFiles <- grep("sample", list.files(file.path(dir_data, 'htseq-count')), 
                    value=TRUE)
metadata <- read.csv(file.path(dir_data, 'metadata.csv'))
sampleTable <- data.frame(sampleName = metadata$sample,
                          fileName = sampleFiles,
                          condition = metadata$group)
sampleTable$condition <- factor(sampleTable$condition)
#---- Gene annotations ---------------------------------------------------------
annotation <- read.csv(file.path(dir_data, 'gencode_GRCh38_v37.csv'))
#---- Build DESeqDataSet -------------------------------------------------------
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = file.path(dir_data, 'htseq-count'),
                                       design= ~ condition)
#---- Raw count ----------------------------------------------------------------
raw_count <- counts(ddsHTSeq) %>%
  as.data.frame() %>%
  mutate(gene_id = row.names(.)) %>%
  inner_join(annotation, by='gene_id') %>%
  select(gene_id, gene_name, gene_type, everything())
write.csv(raw_count, file.path(dir_output, 'raw_count.csv'), row.names = FALSE)
#---- Counts Per Million -------------------------------------------------------
# logcpm
logcpm <- cpm(counts(ddsHTSeq), log=TRUE) %>%
  as.data.frame() %>%
  mutate(gene_id = row.names(.)) %>%
  inner_join(annotation, by='gene_id') %>%
  select(gene_id, gene_name, gene_type, everything())
write.csv(logcpm, file.path(dir_output, 'logcpm.csv'), row.names = FALSE)
# Z scale
z_logcpm <- t(scale(t(cpm(counts(ddsHTSeq), log=TRUE)))) %>%
  as.data.frame() %>%
  mutate(gene_id = row.names(.)) %>%
  inner_join(annotation, by='gene_id') %>%
  select(gene_id, gene_name, gene_type, everything())
write.csv(z_logcpm, file.path(dir_output, 'z_logcpm.csv'))
#---- PCA ----------------------------------------------------------------------
pca_result <- prcomp(t(z_logcpm %>% 
                         select(-c(gene_id, gene_name, gene_type))), 
                     center = TRUE, scale. = FALSE)
pca_data <- data.frame(PC1 = pca_result$x[,1], PC2 = pca_result$x[,2], 
                       Label = sampleTable$condition)
ggplot(pca_data, aes(x = PC1, y = PC2, color = Label)) +
  geom_point() +
  theme_bw() +
  labs(x = "Principal Component 1", 
       y = "Principal Component 2", 
       color = "Group") +
  theme(axis.title = element_text(size = 14), # Axis title
        axis.text = element_text(size = 14), # Axis text
        legend.title = element_text(size = 14), # Legend title
        legend.text = element_text(size = 14) # Legend text
  )
ggsave(file.path(dir_output, 'PCA.png'), width = 8, height = 5)
#---- Differential analysis ----------------------------------------------------
ddsHTSeq <- DESeq(ddsHTSeq)
res <- results(ddsHTSeq, contrast = c('condition', 'group_a', 'group_b'))
res <- as.data.frame(res) %>%
  mutate(gene_id = row.names(res)) %>%
  na.omit() %>%
  mutate(significant = case_when(
    padj < 0.05 & log2FoldChange > 1  ~ 'Up',
    padj < 0.05 & log2FoldChange < -1 ~ 'Down',
    TRUE ~ 'NS'
    )) %>%
  inner_join(annotation, by='gene_id') %>%
  select(gene_id, gene_name, gene_type, everything()) %>%
  arrange(padj)
write.csv(res, file.path(dir_output, 'differential_analysis.csv'))
#---- GSEA ---------------------------------------------------------------------
rnk <- res %>%
  select(gene_name, stat) %>%
  arrange(desc(stat))
write.table(rnk, file.path(dir_output, 'gsea.rnk'), row.names = FALSE, 
            col.names = FALSE, sep = '\t')