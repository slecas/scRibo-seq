#!/usr/bin/env nextflow

process UMI_DEDUP{
    publishDir './scRibo_results', mode: 'copy', pattern:'* _dedup.bam'

    conda 'bioconda::umi_tools bioconda::samtools'
    input:
    path (align_bam)

    output:
    path('*_dedup.bam'), emit:genome_dedup
    tuple val(align_bam.baseName),path('*_dedup.log'), emit: dedup_qc

    publishDir "bam_outs/", mode: 'copy',pattern: '*_dedup.bam'

    script:
    """
    samtools sort -o ${align_bam.baseName}.sorted.bam ${align_bam}

    samtools index ${align_bam.baseName}.sorted.bam

    umi_tools dedup \
    -I ${align_bam.baseName}.sorted.bam\
    -S ${align_bam.baseName}_dedup.bam \
    -L ${align_bam.baseName}_dedup.log \
    --extract-umi-method tag \
    --umi-tag=UR \
    --cell-tag=CB \
    --per-cell \
    --spliced-is-unique \
    --read-length \
    --per-gene \
    --per-contig \
    --no-sort-output
    """
}