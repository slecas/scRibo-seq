#!/usr/bin/env nextflow

// samtools faidx module for scibo pipeline

process SAMTOOLS_FAIDX {
	conda "/home/cashel/projects/nextflow-first/environment.yml"

	input:
	file genome

	output:
	path '*'

	script:
	"""
	samtools faidx ${genome}
	"""
}

