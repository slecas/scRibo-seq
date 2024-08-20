# scriboseqtest
Basic pipeline for processing single cell ribosome profiling (scRibo-seq) data produced by Michael VanInsberghe in 2021 (https://doi.org/10.1038/s41586-021-03887-4). The pipe uses mainly standard bioinformatics tools, with some custom scripts adapted from VanInsberghe's work. 

## Usage
The pipeline is configured to work on fastq files in the format; read 1 <- UMI-read-adpater, read 2 <- Cell barcode (CB). Other inputs required are; a reference genome file, GTF annotation file, CB whitelist, fasta file of sequences to exclude from genome (rRNA etc.). There are also optional inputs that will speed up processing if available; a STAR genome index, a masked genome file and a pickle annotation file.

The pipeline looks for these in the reference directory. To put in any custom files, just direct the relevant input parameter to the desired file path. A small set of test reads are kept in the prac_reads directory. Running the pipeline with test.config from the conf directory will run the pipeline on these reads to make sure it is funcioning.

The pipeline assumes read files containing SRR numbers for grouping. In the case of custom read file names, some small tweaks may be required in the multqc subworkflow.

The pipeline outputs the raw count files from STAR, the deduplicated transcriptome file, a CSV file with random forest predicted p-sites for each read and multiqc reports for each read. These outputs are published in the scRibo-seq directory after the run is complete.

