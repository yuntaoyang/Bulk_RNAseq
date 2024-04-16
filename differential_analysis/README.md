# Differential analysis
This is an example of performing differential analysis using R.
* Check `differential_analysis.R` .
* R requirements: `DESeq2, edgeR, dplyr, ggplot2`.
* Input
    * Output files from htseq-count `./data/htseq-count`.
    * Metadata that includes group information `./data/metadata`.
    * Gene annotation file `gencode_GRCh38_v37.csv`.
* Output
    * PCA plot `PCA.png`.
    * Raw count matrix `raw_count.csv`.
    * Log(Count Per Million) matrix `logcpm.csv`.
    * Z scale Log(Count Per Million) matrix `z_logcpm.csv`.
    * Differentially expressed genes from DESeq2 `differential_analysis.csv`.
    * The input ranked list for GSEA `gsea.rnk`.