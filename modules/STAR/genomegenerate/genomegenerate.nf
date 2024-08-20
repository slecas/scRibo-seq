#!/usr/bin/env nextflow

process STAR_INDEX {
    conda 'bioconda::star=2.7.6a' 

    input:
    path (con_free_genome)
    path (annotation) 

    output:
    path ('gengen'), emit: index


    script:
    //replace everything bar the inputs with params in config

    """
    mkdir gengen

    STAR \
        --runMode genomeGenerate \
        --runThreadN 12 \
        --genomeDir gengen/ \
        --genomeFastaFiles ${con_free_genome} \
        --limitGenomeGenerateRAM 75161927680 \
        --sjdbGTFfile ${annotation} \
        --sjdbOverhang 50 
    """
}

