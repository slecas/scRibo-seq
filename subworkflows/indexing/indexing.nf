#!/usr/bin/env nextflow


include { STAR_INDEX } from '/home/cashel/projects/nextflow-first/modules/STAR/genomegenerate/genomegenerate.nf'
include { CONTAM_REMOVAL } from '/home/cashel/projects/nextflow-first/modules/Bowtie/contamremoval.nf'
include { EXPAND_ANNOTATION } from '/home/cashel/projects/nextflow-first/modules/annotationtools/expandannotation.nf'

workflow INDEXING {

    take:
    genome
    annotations
    contaminants

    main:
    if(params.maskgenome && params.starIndex && params.ex_annotation){
        star_index = Channel.fromPath(params.starIndex)
        p_annot = Channel.fromPath(params.ex_annotation)
        ref_genome = Channel.fromPath(params.maskgenome)
    }

    else if(params.maskgenome && params.starIndex && !params.ex_annotation){
        EXPAND_ANNOTATION(annotations)

        star_index = Channel.fromPath(params.starIndex)
        p_annot = EXPAND_ANNOTATION.out.pickle_annot
        ref_genome = Channel.fromPath(params.maskgenome)
    }
    else if(params.maskgenome && !params.starIndex && params.ex_annotation){
        STAR_INDEX(params.maskgenome, annotations)

        star_index = STAR_INDEX.out.index
        p_annot = Channel.fronPath(params.ex_annotation)
        ref_genome = Channel.fromPath(params.maskgenome)
    }

    else if(!params.maskgenome && params.starIndex && params.ex_annotation){
        CONTAM_REMOVAL(genome,contaminants)

        star_index = Channel.fromPath(params.starIndex)
        p_annot = Channel.fronPath(params.ex_annotation)
        ref_genome = CONTAM_REMOVAL.out.reference_genome
    }

    else if(params.starIndex && !params.maskgenome && !params.ex_annotation){
        CONTAM_REMOVAL(genome,contaminants)
        EXPAND_ANNOTATION(annotations)

        star_index = Channel.fromPath(params.starIndex)
        p_annot = EXPAND_ANNOTATION.out.pickle_annot
        ref_genome = CONTAM_REMOVAL.out.reference_genome
    }

    else if(!params.starIndex && !params.maskgenome && params.ex_annotation){
        CONTAM_REMOVAL(genome,contaminants)
        STAR_INDEX(CONTAM_REMOVAL.out.reference_genome, annotations)
        

        star_index = STAR_INDEX.out.index
        p_annot = Channel.fronPath(params.ex_annotation)
        ref_genome = CONTAM_REMOVAL.out.reference_genome
    }

     else if(!params.starIndex && params.maskgenome && !params.ex_annotation){
        EXPAND_ANNOTATION(annotations)
        STAR_INDEX(params.maskgenome, annotations)
        

        star_index = STAR_INDEX.out.index
        p_annot = EXPAND_ANNOTATION.out.pickle_annot
        ref_genome = Channel.fromPath(params.maskgenome)
    }

    else{
        CONTAM_REMOVAL(genome,contaminants)
        STAR_INDEX(CONTAM_REMOVAL.out.reference_genome, annotations)
        EXPAND_ANNOTATION(annotations)

        star_index = STAR_INDEX.out.index
        p_annot = EXPAND_ANNOTATION.out.pickle_annot
        ref_genome = CONTAM_REMOVAL.out.reference_genome
    }
    emit:
    star_index
    p_annot
    ref_genome 
}