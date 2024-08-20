#!/usr/bin/env nextflow

// module to provide both eukaryotic and mitochondrial trnascan data

process TRNASCAN_SE_MIT {
    publishDir = './'

    conda "bioconda::trnascan-se=2.0.7"

    input: 
    file genome
    path trnatobed

    output:
    path '*_bp.bed' 

    script:
    //run tRNAscan-se
    //replace thread with a param in config
    """
    tRNAscan-SE \\
    -M vert \\
    -Q \\
    -o# \\
    -f# \\
    -m# \\
    -b# \\
    -a# \\
    -l# \\
    --brief \\
    --thread 12 \\
    -p ${genome.baseName}.Mitochondrial_tRNAs \\
    ${genome}

    python3 ${trnatobed} ${genome.baseName}.Mitochondrial_tRNAs.out ${genome.baseName}.Mitochondrial_tRNAs.bed  > ${genome.baseName}.Mitochondrial_tRNAs_bp.bed
    """
}
workflow {
    // establish channel to genome and the python script
    genome = Channel.fromPath('/home/cashel/projects/run_res/parse_reads/sequence.fasta')
    bedscript = Channel.fromPath('/home/cashel/projects/nextflow-first/scripts/trnatobed.py')

    //run process 
    trnascan_ch = TRNASCAN_SE_MIT(genome,bedscript)

    //someway of getting the fucking output
    trnascan_ch.collectFile(name: 'mit_trna_bp.bed', newLine: true, storeDir: 'mit_trna_results').view()

}
