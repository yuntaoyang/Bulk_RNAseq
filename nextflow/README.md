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
* Run the workflow
    * Pull the docker image from Docker Hub `docker pull yuntaoyang1995/rnaseq:yuntaoyang1995/rnaseq:v2`.
    * Revise `nextflow.config` for parameters and computing resources.
    * Run the workflow using `nextflow run rnaseq.nf -profile custom`.