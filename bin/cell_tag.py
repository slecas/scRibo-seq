#!/usr/bin/env python

import pysam
import argparse
import pandas as pd
import os
import tempfile
import shutil

def is_directory(path):
    """Check if the provided path is a valid directory."""
    if os.path.isdir(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"'{path}' is not a valid directory")

def process_file(filepath):
    """Function to process each file."""
    print(f"Processing file: {filepath}")
    # Add your file processing code here

def add_cycle_tag(input_bam, identity, temp_bam):
    """Add a cell cycle tag to reads in the BAM file and write to a temporary BAM file."""
    with pysam.AlignmentFile(input_bam, "rb") as bam_in:
        bam_samp = str(bam_in.filename).split('/')[-1][9:11] # Extract sample name from filename
        # Read and filter the identity file for the specific sample
        id_file = pd.read_csv(identity, sep='\t')
        id_sub = id_file[id_file['sample'] == int(bam_samp)]
        barcode_dict = id_sub.set_index('barcode')['phase'].to_dict()
        with pysam.AlignmentFile(temp_bam, "wb", header=bam_in.header) as bam_out:
            for read in bam_in:
                barcode = read.get_tag('CB') if read.has_tag('CB') else None

                if barcode and barcode in barcode_dict.keys():
                    read.set_tag('SC', barcode_dict[barcode])
                bam_out.write(read)

def concatenate_bams(bam_files, output_bam):
    """Concatenate all BAM files into a single output BAM file."""
    with pysam.AlignmentFile(output_bam, "wb", header=pysam.AlignmentFile(bam_files[0], "rb").header) as bam_out:
        for bam_file in bam_files:
            with pysam.AlignmentFile(bam_file, "rb") as bam_in:
                for read in bam_in:
                    bam_out.write(read)

def main():
    parser = argparse.ArgumentParser(description="Process to add a cell cycle tag to reads in a bam file")
    parser.add_argument('-ct','--cellstage', type=str, help="Path to tsv structured as sample:barcode:cycle stage")
    parser.add_argument('-d','--in_dir', type=is_directory, help="Path to the directory BAM file without CB tags.")
    parser.add_argument('-o','--output', type=str, help="Path to the output BAM file with added CB tags.")
    
    args = parser.parse_args()
    
    input_dir = args.in_dir
    output_bam_path = args.output if args.output else os.path.join(os.getcwd(), 'concatenated_output.bam')
    
    # Get the list of BAM files in the directory
    bam_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith('.bam')]
    
        # Temporary files for each processed BAM
    temp_bams = []
    
    for bam_file in bam_files:
        # Create a temporary file
        temp_fd, temp_bam_path = tempfile.mkstemp(suffix='.bam')
        os.close(temp_fd)
        temp_bams.append(temp_bam_path)
        
        print(f"Processing file: {bam_file}")
        add_cycle_tag(bam_file, args.cellstage, temp_bam_path)
    
    # Concatenate all temporary BAM files into the final output BAM file
    print(f"Concatenating BAM files into {output_bam_path}")
    concatenate_bams(temp_bams, output_bam_path)
    
    # Cleanup temporary files
    for temp_bam_path in temp_bams:
        os.remove(temp_bam_path)

    print(f"All BAM files processed and concatenated into {output_bam_path}")

if __name__ == "__main__":
    main() 
