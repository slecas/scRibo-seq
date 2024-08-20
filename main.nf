#!/usr/bin/env nextflow

// ok i mean might aswell start with the pre-processing part
// so this file needs to effetively just be the workflow, with all processes imported

println """\
  scRIBO-seq processing pipeline
  ==============================
  reads       :${params.reads}
  cblen       :${params.cblen}
  umilen      :${params.umilen}
  annotations :${params.annotations}
  contaminants:${params.contamination}
  genome      :${params.genome}
  whitelist   :${params.whitelist}
  """
  .stripIndent()

// including subworkflows
include { ALIGNMENT } from '/home/cashel/projects/nextflow-first/subworkflows/alignment/alignment.nf'
include { PREPROCESS } from '/home/cashel/projects/nextflow-first/subworkflows/preprocess/preprocess.nf'
include { INDEXING } from '/home/cashel/projects/nextflow-first/subworkflows/indexing/indexing.nf'
include { MULTIQC } from '/home/cashel/projects/nextflow-first/subworkflows/multqc/multqc.nf'

// establishing channels for parameters
gen_ch = Channel.fromPath(params.genome)
an_ch = Channel.fromPath(params.annotations)
con_ch = Channel.fromPath(params.contamination)
read_ch = Channel.fromFilePairs(params.reads)
model = Channel.fromPath(params.trained_model)
whitelist = Channel.fromPath(params.whitelist)

workflow SCRIBOSEQ {
    main:
    PREPROCESS(read_ch,params.cblen,params.umilen)
    INDEXING (gen_ch,an_ch,con_ch)

    ref = INDEXING.out.ref_genome.collect()
    starIND =INDEXING.out.star_index.collect()
    p_annot = INDEXING.out.p_annot.collect()

    ALIGNMENT(starIND,whitelist,params.annotations,p_annot,ref,model,PREPROCESS.out.trimmed_reads)
    
    MULTIQC(PREPROCESS.out.raw_fastqc,PREPROCESS.out.trim_fastqc,PREPROCESS.out.cut_qc,ALIGNMENT.out.star_qc,ALIGNMENT.out.dedup_qc)

  emit:
  ALIGNMENT.out.matrix_results
  ALIGNMENT.out.p_site_dict
  MULTIQC.out.m_html
}

