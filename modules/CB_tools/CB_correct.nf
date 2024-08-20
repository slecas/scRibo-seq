#!/usr/bin/env nextflow

process CB_CORRECT {
    conda'bioconda::pysam conda-forge::pandas=1.2.2'

    input:
    path(bam_trans)
    path(bam_gen)

    output:
    path('*_cb_corr.bam'), emit: corr_bam

    script:
    """
    cb_lookup.py ${bam_gen} ${bam_trans} ${bam_trans.baseName}_cb_corr.bam
    """

}