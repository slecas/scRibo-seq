#!/usr/bin/env nextflow

process EXPAND_ANNOTATION {
    conda'conda-forge::python conda-forge:pandas conda-forge::feather-format bioconda::gffutils'

    input:
    path (gtf_file)

    output:
    path ('*.annotations.pickle'), emit: pickle_annot

    script:
    """
    expandannot.py -g ${gtf_file}
    """
}
