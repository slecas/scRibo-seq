#!/usr/bin/env nextflow

process UMI_DEDUP{
    conda 'bioconda::umi_tools bioconda::samtools'
    input:
    path (align_bam)

    output:
    path ('*_dedup.bam'), emit:genome_dedup
    path('*_dedup.log'), emit: dedup_qc

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

    //change to CR for transcriptome
    //samtools view -h ${align_bam} | awk 'BEGIN {OFS="\\t"} /^@/ {print; next} {
      //cb_empty = 1;
      //ub_empty = 1;
      //for (i = 12; i <= NF; i++) {
        //if (\$i ~ /^CR:Z:/) cb_empty = (\$i == "CR:Z:-");
        //if (\$i ~ /^UR:Z:/) ub_empty = (\$i == "UR:Z:-");
      //}
      //if (!cb_empty && !ub_empty) print
    //}' | samtools view -bS - > ${align_bam.baseName}_valid.bam