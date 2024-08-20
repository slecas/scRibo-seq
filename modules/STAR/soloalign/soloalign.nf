#!/usr/bin/env nextflow
process SOLO_ALIGN {
    conda 'bioconda::STAR'
    maxForks 1

    input:
    tuple path(star_index),val(name),path(read1),path(read2),path(whitelist)

    output:
    path ('*.sortedByCoord.out.bam'), emit:aligned_genome
    path ('*.toTranscriptome*'), emit:aligned_transcriptome
    path ('*Solo.out'),emit:matrix_results

    publishDir ".", mode: 'copy',pattern: '*Solo.out'

    script:
    """

    STAR \
    --outFileNamePrefix ${name} \
    --genomeDir ${star_index} \
    --runThreadN 10 \
    --readFilesIn ${read1} ${read2}  \
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
    --outSAMattributes NH HI AS nM NM MD jM jI MC ch CR CB UR UY GX GN CY UY sS sQ sM
    """
}


