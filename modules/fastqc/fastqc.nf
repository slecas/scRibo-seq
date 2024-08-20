#!/usr/bin/env nextflow 

process RAWFASTQC{
    conda 'bioconda::fastqc'

    input:
    tuple val(name), path (reads)

    output:
    path ('*'), emit: fastqc

    script:
    """
    fastqc -t 12 -f fastq -q ${reads[0]} --nogroup
    """
}
