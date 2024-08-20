#!/usr/bin/env nextlfow
process CUTADAPT {
    conda 'bioconda::cutadapt conda-forge::python=3.8.2'

    input:
    tuple val(name), path (reads)

    output:
    tuple val(name),path("*.trimmed.fastq"), emit: trimmed_reads
    path ("*.trim_report.txt"), emit: multiqc_report

    script:
    """
    cutadapt \
        --cores=12 \
        -m 15 \
        -a TGGAATTCTCGGGT \
        -o ${reads[0].baseName}_R1.trimmed.fastq \
        -p ${reads[1].baseName}_R2.trimmed.fastq \
        ${reads[0]} ${reads[1]} > ${reads[0].baseName}.trim_report.txt
    """
}

//mv ${reads[1]} ${reads[1].baseName}_R2.trimmed.fastq