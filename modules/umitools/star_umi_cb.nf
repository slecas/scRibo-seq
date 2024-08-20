#!/usr/bin/env nextflow

process STAR_UMI_CB{
    conda'conda-forge::python bioconda::pysam'

    input:
    tuple val(name),path (reads)
    val (umi_length)

    output:
    tuple val(name),path('*_cbumi.fastq'), emit: cb_umi_reads

    script:
    """
    star_umi_cb.py ${reads[0]} ${reads[1]} ${reads[0].baseName}_cbumi.fastq ${reads[1].baseName}_cbumi.fastq 
    """
}
