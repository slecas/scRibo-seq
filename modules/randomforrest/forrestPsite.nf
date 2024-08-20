#!/usr/bin/env nextflow 

process RANDOM_FORREST {
    conda'conda-forge::r-base r-dplyr r-tidyr r-data.table r-optparse r-feather r-forcats r-ranger lapack r-r.utils conda-forge::gzip'

    input:
    tuple path (bam_file),path (model)

    output:
    path ('*.csv.gz'), emit: psite_adjusted_bam
    path ('*_dict.csv.gz'), emit: psite_dictionary

    script:
    """
    applymodel.R \
    -m ${model} \
    -r ${bam_file}

    cut -f 10,47 > ${bam_file.baseName}_dict.csv | gzip
    """ 
}