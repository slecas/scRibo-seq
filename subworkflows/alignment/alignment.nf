#!/usr/bin/env nextflow

//params.star_index = ''
//params.trim_reads = ''
//params.whitelist = ''
//params.annotation = ''

include { SOLO_ALIGN } from '/home/cashel/projects/nextflow-first/modules/STAR/soloalign/soloalign.nf'
include { UMI_DEDUP } from '/home/cashel/projects/nextflow-first/modules/umitools/umi_dedup.nf'
include { CB_CORRECT } from '/home/cashel/projects/nextflow-first/modules/CB_tools/CB_correct.nf'
include { PARSE_READS } from '/home/cashel/projects/nextflow-first/modules/annotationtools/parsereads.nf'
include { RANDOM_FORREST } from '/home/cashel/projects/nextflow-first/modules/randomforrest/forrestPsite.nf'

workflow ALIGNMENT {
    take: 
    star_index
    whitelist
    annotation
    pick_annot
    genome
    forrest_model
    trim_reads

    main:
    index = star_index.collect()
    white = whitelist.collect()
    SOLO_ALIGN(index,trim_reads,white)

    CB_CORRECT(SOLO_ALIGN.out.aligned_transcriptome,SOLO_ALIGN.out.aligned_genome)
    UMI_DEDUP(CB_CORRECT.out.corr_bam)

    p_anno = pick_annot.collect()
    gen = genome.collect()
    PARSE_READS(UMI_DEDUP.out.genome_dedup,p_anno,gen)

    model = forrest_model.collect()
    RANDOM_FORREST(PARSE_READS.out.csv_out,model)

    emit:
    matrix_results = SOLO_ALIGN.out.matrix_results
    star_qc = SOLO_ALIGN.out.star_qc
    dedup_qc = UMI_DEDUP.out.dedup_qc
    p_site_dict = RANDOM_FORREST.out.psite_dictionary
}

