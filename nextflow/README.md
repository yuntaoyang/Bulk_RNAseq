# Nextflow
This is a workflow that can automatically perform adapter trimming, quality control, and alignment.
* Check `rnaseq.nf` .
* Requirements: `nextflow, trim_galore, STAR, multiqc`.
* Input data
    * Metadata `./data/metadata.csv`.
    * Reference genome `./data/reference_genome.fasta.gz`.
    * Gene annotation `./data/gene_annotation.gtf.gz`.
* Output data
    * Raw reads count `./output/sample_id_ReadsPerGene.out.tab`.
    * MultiQC report `./output/multiqc_report.html`.