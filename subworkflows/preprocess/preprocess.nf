#!/usr/bin/env nextflow

include { CUTADAPT } from '/home/cashel/projects/nextflow-first/modules/fastqc/cutadapt.nf'
include { RAWFASTQC } from '/home/cashel/projects/nextflow-first/modules/fastqc/fastqc.nf'
include { STAR_UMI_CB } from '/home/cashel/projects/nextflow-first/modules/umitools/star_umi_cb.nf'
include { RAWFASTQC as TRIMMEDFASTQC} from '/home/cashel/projects/nextflow-first/modules/fastqc/fastqc.nf'

workflow PREPROCESS{

    take:
    reads
    umilen
    cblen

    main:
    STAR_UMI_CB(reads,umilen)
    CUTADAPT(STAR_UMI_CB.out.cb_umi_reads)
    RAWFASTQC (reads)
    TRIMMEDFASTQC(CUTADAPT.out.trimmed_reads)

    emit:
    trimmed_reads = CUTADAPT.out.trimmed_reads
    cut_qc = CUTADAPT.out.multiqc_report
    raw_fastqc = RAWFASTQC.out.fastqc
    trim_fastqc = TRIMMEDFASTQC.out.fastqc
}