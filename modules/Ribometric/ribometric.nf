#!/usr/bin/env nextflow

process RIBOMETRIC{
    conda'conda-forge:python=3.8'

    input:
    tuple path(annotations),path (bam_file)

    output:
    path ('*.html'), emit:ribo_pipe

    publishDir ".", mode: 'copy',pattern: '*.html'

    script:
    """
    pip install RiboMetric
    
    RiboMetric prepare -g ${annotations}

    RiboMetric run -b ${bam_file} -a ${annotations}.tsv
    """
}