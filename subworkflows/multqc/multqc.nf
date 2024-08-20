#!/usr/bin/env nextflow

include { MULTIQC_END } from '/home/cashel/projects/nextflow-first/modules/multqc/multqc_mod.nf'

workflow MULTIQC {
    take:
    raw_qc
    trim_qc
    cut_qc
    star_qc
    umi_qc

    main:
    raw_q = raw_qc.map { 
        def id = it[0].toString().split('/')[-1].substring(0,11)
        def pathString2 = it[1]  // Convert UnixPath to String
        def pathString1 = it[0]  // Assuming the ID is in the first element of the tuple
        [id,pathString1,pathString2]
    }

    trim_q = trim_qc.map { 
        def id = it[0].toString().split('/')[-1].substring(0,11)
        def pathString2 = it[1]  // Convert UnixPath to String
        def pathString1 = it[0]  // Assuming the ID is in the first element of the tuple
        [id,pathString1,pathString2]
    }

    // Process cut channel
    cut_q = cut_qc.map { 
        def id = it.toString().split('/')[-1].substring(0,11)
        def pathString = it.toString()  // Convert UnixPath to String
        [id,pathString]
    }

    star_q = star_qc.map { 
        def id = it[0].toString().substring(0,11)
        def pathLog = it[1].toString()  // Convert UnixPath to String
        def pathReads = it[2].toString()  // Convert UnixPath to String
        [id, pathLog, pathReads]
    }.view()

    umi_q = umi_qc.map { 
        def pathString = it[1].toString() // Convert UnixPath to String
        def id = it[0].toString().substring(0,11)
        [id, pathString]
    }.view()

    m_input = raw_q.join(trim_q, by : 0).join(cut_q,by:0).join(star_q, by : 0).join(umi_q, by : 0)

    MULTIQC_END(m_input)

    emit:
    m_html = MULTIQC_END.out.html

}