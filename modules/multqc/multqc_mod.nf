#!/usr/bin/env nextflow

process MULTIQC_END{
    publishDir './scRibo_results',mode: 'copy',pattern:'*.html'
    conda'bioconda::multiqc'

    input:
    tuple val(name),path(raw1),path(raw2),path(trim1),path(trim2),path(cut),path(star1),path(star2),path(umi)

    output:
    path("multiqc_report.html"),emit:html

    script:
    """
    multiqc .
    """
}