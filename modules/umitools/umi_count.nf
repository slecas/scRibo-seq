#!/usr/bin/env nextflow

//this is actually feature counts but hey

include { SAMBAMBA_SORT } from '/home/cashel/projects/nextflow-first/modules/samtools/samtools_sort.nf'

process FEATURE_COUNTS {
    conda'bioconda::subread'

    input:
    path (sorted_bam)
    path (annotation)

    output:
    path ('*.featureCounts.bam'), emit:count_bam

    script:
    """
    featureCounts -a ${annotation} -o counts.txt -R BAM ${sorted_bam}
    """
}
