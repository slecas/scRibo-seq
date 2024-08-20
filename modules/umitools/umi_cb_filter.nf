#!/usr/bin/env nextflow

process CB_FILTER {
    conda'bioconda::umis'

    input:
    path(bam_file)
    path(whitelist)

    output:
    path('*_CB.bam'), emit:filtered_transcriptome

    script:
    """
    umis cb_filter \
    --stdin=${bam_file} \
    --stdout=${bam_file.baseName}_CB.bam \
    --whitelist=${whitelist} \
    --nedit 1
    """
}
	
