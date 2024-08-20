#!/usr/bin/env nextflow

process PARSE_READS {
    conda'conda-forge::python conda-forge::pandas bioconda::gffutils bioconda::pysam conda-forge::feather-format bioconda::samtools'

    input:
    tuple path (bam_file),path (pick_annot),path(genome)

    output:                                                                                                                                                                                                                                                                                                                                        
    path('*.csv.gz'), emit:csv_out

    script:
    """
    samtools sort ${bam_file} -o ${bam_file}.sorted.bam
    samtools index ${bam_file}.sorted.bam
    
    parsereads.py -t ${pick_annot} -f ${genome} -b ${bam_file}.sorted.bam -ba ${bam_file}.sorted.bam.bai --id
    """
}