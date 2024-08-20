#!/usr/bin/env nextflow 

process RANDOM_FORREST {
    publishDir './scRibo_results',mode: 'copy',pattern:'*_dict.csv'
    // a not: this model is VERY memory intensive. its limited to one fork here
    // and the R script applying it has been set into chunks. If issue arises
    // changing the chunk size in the script is yout best bet


    conda'conda-forge::r-base r-dplyr r-tidyr r-data.table r-optparse r-feather r-forcats r-ranger lapack r-r.utils'
    maxForks 1

    input:
    path (bam_file)
    path (model)

    output:
    path ('*.csv'), emit: psite_adjusted_bam
    path ('*_dict.csv'), emit: psite_dictionary

    script:
    """
    applymodel.R \
    -m ${model} \
    -r ${bam_file} \

    cut -d ',' -f 11,47 ${bam_file.baseName}.predicted.csv > ${bam_file.baseName}_dict.csv
    """ 
}