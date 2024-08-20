#!/usr/bin/env nextflow

process CREATE_GEN{

    conda 'bioconda::bedtools bioconda::samtools'


    publishDir = './genomemask/'

    input:
    file euk_mask
    file mit_mask
    file genome
    file annotation
    file genome_index
    path expandbed
    path appendfasta
    path collapseSequences
    path mergeGTF

    output: 
    path ('*')

    script:
    """
    ### tRNAs
    #mask all tRNAs in genome
    cat ${euk_mask} ${mit_mask} | sort -k 1,1 -k2,2n > ${genome.baseName}.tRNA_masked_regions.bed
    bedtools maskfasta -fi ${genome} -fo ${genome.baseName}.tRNA_masked.fa -bed ${genome.baseName}.tRNA_masked_regions.bed

    #pre-tRNAs
    grep -v "pseudo\\s" ${euk_mask} | python3 ${expandbed} --slop 50 --index ${genome_index} > ${euk_mask.baseName}.pre-tRNAs.bed

    #mit-tRNAs
    sed '/^MT\\|^chrM/!d' ${mit_mask} > ${mit_mask.baseName}.only_MT.bed
    grep -v "pseudo\\s" ${mit_mask.baseName}.only_MT.bed | python3 ${expandbed} --slop 50 --index ${genome_index} > ${mit_mask.baseName}.pre-tRNAs.bed

    #combine both
    cat ${euk_mask.baseName}.pre-tRNAs.bed ${mit_mask.baseName}.pre-tRNAs.bed | sort -k 1,1 -k2,2n > ${genome.baseName}.pre-tRNAs.bed

    #extract the pre-tRNA sequences
    bedtools getfasta -name -split -s -fi ${genome} -bed ${genome.baseName}.pre-tRNAs.bed -fo ${genome.baseName}.pre-tRNAs.fa

    #mature tRNAs
    cat ${euk_mask} ${mit_mask.baseName}.only_MT.bed | sort -k 1,1 -k2,2n > ${genome.baseName}.mature-tRNAs.bed
    bedtools getfasta -name -split -s -fi ${genome} -bed ${genome.baseName}.mature-tRNAs.bed | python3 ${appendfasta} --append cca > ${genome.baseName}.mature-tRNAs.fa

    #cluster stuff
    ## cluster mature tRNAs
    python3 ${collapseSequences} ${genome.baseName}.mature-tRNAs.fa > ${genome.baseName}.mature-tRNAs.cluster.fa
    # produces cluster_info.gtf
    mv cluster_info.gtf ${genome.baseName}.tRNACluster.gtf
    samtools faidx ${genome.baseName}.mature-tRNAs.cluster.fa

    #final genome assembly
    cat ${genome.baseName}.tRNA_masked.fa ${genome.baseName}.mature-tRNAs.cluster.fa > ${genome.baseName}.tRNA-masked-withmature.fa
    samtools faidx ${genome.baseName}.tRNA-masked-withmature.fa

    #### Annotations
    python3 ${mergeGTF} ${annotation} ${genome.baseName}.tRNACluster.gtf > ${annotation.baseName}.tRNA.gtf
    """
}

workflow{
    euk_mask = Channel.fromPath('/home/cashel/projects/nextflow-first/modules/trnascan/results/trna_bp.bed')
    mit_mask = Channel.fromPath('/home/cashel/projects/nextflow-first/modules/trnascan/mit_trna_results/mit_trna_bp.bed')
    genome = Channel.fromPath('/home/cashel/projects/run_res/parse_reads/sequence.fasta')
    genome_index = Channel.fromPath('/home/cashel/projects/run_res/parse_reads/sequence.fasta.fai')
    expandbed = Channel.fromPath('/home/cashel/projects/nextflow-first/scripts/expandbed.py')
    annotation = Channel.fromPath('/home/cashel/projects/run_res/chr1annot.gtf')
    appendfasta = Channel.fromPath('/home/cashel/projects/nextflow-first/scripts/appendfasta.py')
    collapseSequences = Channel.fromPath('/home/cashel/projects/nextflow-first/scripts/collapseseq.py')
    mergeGTF = Channel.fromPath('/home/cashel/projects/nextflow-first/scripts/mergegtf.py')

    CREATE_GEN(euk_mask,mit_mask,genome,annotation,genome_index,expandbed,appendfasta,collapseSequences,mergeGTF)

    CREATE_GEN.out
}