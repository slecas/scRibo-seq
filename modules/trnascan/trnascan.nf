#!/usr/bin/env nextflow

// module to provide both eukaryotic and mitochondrial trnascan data

process TRNASCAN_SE {

    conda "bioconda::trnascan-se=2.0.7"

    input: 
    file genome
    path bedscript

    output:
    path '*_bp.bed' 

    script:
    //run tRNAscan-se
    // replace thread with a param in config
    """
    tRNAscan-SE \\
    -HQ \\
    -o# \\
    -f# \\
    -s# \\
    -m# \\
    -b# \\
    -a# \\
    -l# \\
    --brief \\
    --thread 12 \\
    -p ${genome.baseName}.Eukaryotic_tRNAs \\
    ${genome}

    python3 ${bedscript} ${genome.baseName}.Eukaryotic_tRNAs.out ${genome.baseName}.Eukaryotic_tRNAs.bed > ${genome.baseName}.Eukaryotic_tRNAs_bp.bed |
    sed '/^MT\\|^chrM/d' ${genome.baseName}.Eukaryotic_tRNAs_bp.bed > ${genome.baseName}.Eukaryotic_tRNAs_bp.bed
    """
}

workflow {
    // establish channel to genome and the python script
    genome = Channel.fromPath('/home/cashel/projects/run_res/parse_reads/sequence.fasta')
    bedscript = Channel.fromPath('/home/cashel/projects/nextflow-first/scripts/trnatobed.py')

    //run process 
    trnascan_ch = TRNASCAN_SE(genome,bedscript)

    //someway of getting the fucking output
    trnascan_ch.collectFile(name: 'euk_trna_bp.bed', newLine: true, storeDir: 'euk_trna_results').view()
}

