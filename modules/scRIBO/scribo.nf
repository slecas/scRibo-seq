#!/usr/bin/ nextflow

// scribo annotation time

process SCRIBO_ANNOTATION{
    conda 'conda-forge::python=3.8.2 conda-forge::numpy=1.20.0 conda-forge::pandas=1.2.2 conda-forge::feather-format=0.4.1 bioconda::gffutils bioconda::pysam'
    publishDir = './'

    input:
    file genome
    file annotation
    path exanno
    path excod

    output:
    path ('*')

    script:
    """
    expandannot.py -g ${annotation}
    expandcodon.py -f ${genome} -t ${annotation.baseName}.annotations.pickle
    """
}

workflow{
    genome = Channel.fromPath('/home/cashel/projects/run_res/parse_reads/sequence.fasta')
    annotation = Channel.fromPath('/home/cashel/projects/run_res/chr1annot.gtf')
    exanno = Channel.fromPath('/home/cashel/projects/nextflow-first/scripts/expandannot.py')
    excod = Channel.fromPath('/home/cashel/projects/nextflow-first/scripts/expandcodon.py')

    SCRIBO_ANNOTATION(genome,annotation,exanno,excod)

    SCRIBO_ANNOTATION.out
}

