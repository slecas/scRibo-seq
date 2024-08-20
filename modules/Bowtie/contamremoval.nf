#!/usr/bin/env nextflow

process CONTAM_REMOVAL {
    conda 'bioconda::bowtie2 bioconda::samtools bioconda::bedtools'

    input:
    path (genome)
    path (contaminants)

    output:
    path ('*_clean.fasta'), emit: reference_genome

    script:
    """
    #build genome index; multiple bt2 files
    bowtie2-build ${genome} ${genome.baseName}_index

    #use these files to create contaminant map
    bowtie2 -x ${genome.baseName}_index -f ${contaminants} -S contaminants.sam

    #samtools manipulation to extract areas in refernce that match contaminant regions
    samtools view -bS contaminants.sam > contaminants.bam
    samtools sort contaminants.bam -o contaminants_sorted.bam
    samtools index contaminants_sorted.bam

    #convert to bed file for easier use
    bedtools bamtobed -i contaminants_sorted.bam > contaminants.bed

    #mask original reference
    bedtools maskfasta -fi ${genome} -bed contaminants.bed -fo ${genome.baseName}_clean.fasta
    """ 
}
