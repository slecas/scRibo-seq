#!/usr/bin/env nextflow 

process SAMBAMBA_SORT{
    conda'bioconda::sambamba'

    input:
    path (bam_file)

    output:
    path ('*.sorted.bam'), emit: sorted_bam

    script:
    """
    sambamba sort \
    -m 2G \
    -t 12 \
    ${bam_file}
    """
}


