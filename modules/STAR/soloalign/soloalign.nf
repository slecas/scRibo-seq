#!/usr/bin/env nextflow

process SOLO_ALIGN {
    publishDir './scRibo_results', mode: 'copy', pattern: '*Solo.out'

    conda 'bioconda::STAR'
    maxForks 1

    input:
    path(star_index)
    tuple val(name),path(reads)
    path(whitelist)

    output:
    tuple val(name), path('*.final.out'),path('*ReadsPerGene.out.tab'), emit: star_qc
    path ('*.sortedByCoord.out.bam'), emit:aligned_genome
    path ('*.toTranscriptome*'), emit:aligned_transcriptome
    path ('*Solo.out'),emit:matrix_results


    script:
    """
    STAR \
    --outFileNamePrefix ${name} \
    --genomeDir ${star_index} \
    --runThreadN 10 \
    --readFilesIn ${reads[0]} ${reads[1]}  \
    --outFilterMultimapNmax 1 \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM GeneCounts \
    --soloType CB_UMI_Simple \
    --soloCBstart 11 \
    --soloCBlen 10 \
    --soloUMIstart 1 \
    --soloUMIlen 10 \
    --soloBarcodeReadLength 0 \
    --soloCBwhitelist ${whitelist} \
    --soloStrand Forward \
    --soloCBmatchWLtype 1MM \
    --soloFeatures Gene GeneFull Velocyto \
    --soloCellFilter None \
    --soloUMIdedup 1MM_Directional \
    --outSAMattributes ${params.star_out}
    """
}


