#!/usr/bin/env nextflow

// Enable DSL2
nextflow.enable.dsl = 2

// Parameters
params.input = "samplesheet.csv"
params.outdir = "results"

// FastQC on raw reads
process fastqc_raw {
    module 'fastqc/0.12.1'
    publishDir "${params.outdir}/fastqc_raw", mode: 'copy'
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.{zip,html}"

    script:
    """
    echo "Running FastQC on raw reads: ${sample_id}"
    fastqc ${reads}
    """
}

// Trimmomatic for quality trimming
process trimmomatic {
    module 'trimmomatic/0.39'
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_*_paired.fastq.gz")
    path "${sample_id}_*_unpaired.fastq.gz"

    script:
    """
    echo "Running Trimmomatic on ${sample_id}"

    trimmomatic PE -threads 2 \\
        ${reads[0]} ${reads[1]} \\
        ${sample_id}_R1_paired.fastq.gz ${sample_id}_R1_unpaired.fastq.gz \\
        ${sample_id}_R2_paired.fastq.gz ${sample_id}_R2_unpaired.fastq.gz \\
        LEADING:3 TRAILING:3 \\
        SLIDINGWINDOW:4:15 MINLEN:36
    """
}

// FastQC on trimmed reads
process fastqc_trimmed {
    module 'fastqc/0.12.1'
    publishDir "${params.outdir}/fastqc_trimmed", mode: 'copy'
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.{zip,html}"

    script:
    """
    echo "Running FastQC on trimmed reads: ${sample_id}"
    fastqc ${reads}
    """
}

// SPAdes genome assembly
process spades_assembly {
    module 'spades/4.2.0'
    publishDir "${params.outdir}/assemblies", mode: 'copy'
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_assembly/contigs.fasta")
    path "${sample_id}_assembly/"

    script:
    """
    echo "Running SPAdes assembly on ${sample_id}"

    spades.py \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -o ${sample_id}_assembly \\
        --threads 2 \\
        --memory 8
    """
}

// Prokka genome annotation
process prokka_annotation {
    module 'prokka/1.14.6'
    publishDir "${params.outdir}/annotation", mode: 'copy'
    tag "$sample_id"

    input:
    tuple val(sample_id), path(contigs)

    output:
    path "${sample_id}_annotation/"

    script:
    """
    echo "Running Prokka annotation on ${sample_id}"

    prokka \\
        --outdir ${sample_id}_annotation \\
        --prefix ${sample_id} \\
        --cpus 2 \\
        --genus Mycobacterium \\
        --species tuberculosis \\
        --kingdom Bacteria \\
        ${contigs}
    """
}

// Main workflow
workflow {
    // Read sample sheet and create channel
    read_pairs_ch = Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def sample = row.sample
            def fastq1 = file(row.fastq_1)
            def fastq2 = file(row.fastq_2)
            return [sample, [fastq1, fastq2]]
        }

    // Run FastQC on raw reads
    fastqc_raw_results = fastqc_raw(read_pairs_ch)
    fastqc_raw_results.view { "Raw FastQC: $it" }

    // Run Trimmomatic for quality trimming
    (trimmed_paired, trimmed_unpaired) = trimmomatic(read_pairs_ch)
    trimmed_paired.view { "Trimmed paired reads: $it" }

    // Run FastQC on trimmed reads
    fastqc_trimmed_results = fastqc_trimmed(trimmed_paired)
    fastqc_trimmed_results.view { "Trimmed FastQC: $it" }

    // Run SPAdes assembly
    (assembly_contigs, assembly_dir) = spades_assembly(trimmed_paired)
    assembly_contigs.view { "Assembly contigs: $it" }

    // Run Prokka annotation
    annotations = prokka_annotation(assembly_contigs)
    annotations.view { "Annotation: $it" }
}
