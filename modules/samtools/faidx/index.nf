#!/usr/bin/env nextflow

// samtools index module for scibo pipeline

process SAMTOOLS_INDEX {
	conda "bioconda::samtools"

	input:
	path (raw_file)

	output:
	path ('*.bai'), emit:index_file

	script:
	"""
	samtools index ${raw_file}
	"""
}

