#!/usr/bin/env nextflow

process UMITOOLS_COUNT{
    conda'bioconda::umi_tools'

    input:
    path (count_bam)
    path (index_file)

    output:
    path ('*.tsv'),emit:count_res

    script:
    """
    umi_tools count \
    -I ${count_bam} \
    -S scribo_count.tsv \
    --extract-umi-method tag \
    --umi-tag=UR \
    --cell-tag=CR \
    --per-gene \
    --gene-tag=XT \
    --assigned-status-tag=XS \
    """
}