#!/usr/bin/env nextflow
// Parameters
params.metadata = "$projectDir/data/metadata.csv"
params.genomeFile = "$projectDir/data/reference_genome.fasta.gz"
params.gtfFile = "$projectDir/data/gene_annotation.gtf.gz"
params.outdir = "$projectDir/output/"

// TrimGalore process to remove adapter and quality control
process TrimGalore {
    tag "TrimGalore on $sample_id"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), val(group), 
          file(read1), file(read2)

    output:
    tuple val(sample_id),
          path("${sample_id}_val_1.fq.gz"), 
          path("${sample_id}_val_2.fq.gz"),
          path("${sample_id}_val_1_fastqc.zip"), 
          path("${sample_id}_val_2_fastqc.zip")

    script:
    """
    trim_galore --fastqc --gzip --paired --basename ${sample_id} ${read1} ${read2}
    """
}

// STAR index
process StarIndex {
    tag "STAR index"

    input:
    path genome
    path gtf

    output:
    path "star_index"

    script:
    """
    zcat $genome > reference_genome.fasta
    zcat $gtf > gene_annotation.gtf
    STAR --runThreadN 16 \
         --runMode genomeGenerate \
         --genomeDir star_index \
         --genomeFastaFiles reference_genome.fasta \
         --sjdbGTFfile gene_annotation.gtf
    """
}

// STAR alignment
process StarAlign {
    tag "STAR on $sample_id"
    publishDir params.outdir, mode:'copy'

    input:
    path star_index
    tuple val(sample_id),
          path(read1), path(read2)

    output:
    tuple path("${sample_id}_ReadsPerGene.out.tab"),
          path("${sample_id}_Log.final.out")

    script:
    """
    STAR --genomeDir $star_index \
         --readFilesIn ${read1} ${read2} \
         --runThreadN 16 \
         --outSAMtype BAM SortedByCoordinate \
         --readFilesCommand zcat \
         --quantMode GeneCounts \
         --outFileNamePrefix ${sample_id}_
    """
}

// MultiQC process to aggregate results from Trim-Galore and STAR
process MultiQC {
    tag "MultiQC report"
    publishDir params.outdir, mode:'copy'

    input:
    path qc_files

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc .
    """
}

workflow {
    // Read CSV and get metadata
    Channel
        .fromPath(params.metadata)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row -> tuple(row.sample_id, row.group, file(row.file_1), file(row.file_2)) }
        .set { samples_ch }
    // Trim Galore
    TrimGalore(samples_ch)
        .set { trimmed_ch }
    // STAR
    trimmed_ch
        .map { sample_id, read1, read2, fastqc_zip1, fastqc_zip2 ->
            tuple(sample_id, read1, read2)
        }
        .set { trimmed_ch_star }
    star_index_ch = StarIndex(params.genomeFile, params.gtfFile)
    align_ch = StarAlign(star_index_ch, trimmed_ch_star)
    // MultiQC
    trimmed_ch
        .map { sample_id, read1, read2, fastqc_zip1, fastqc_zip2 ->
            tuple(fastqc_zip1, fastqc_zip2)
        }
        .set { fastqc_files }
    align_ch
        .map { reads_per_gene, log_final_out -> 
            tuple(log_final_out)
        }
        .set { star_log_files }
    MultiQC(fastqc_files.mix(star_log_files).collect())
}