#!/usr/bin/env python

import pysam
import argparse
import pandas as pd
import os
import tempfile
import shutil

def get_unique_sc_tags(input_bam):
    """Get a set of unique SC tags from the BAM file."""
    sc_tags = set()
    with pysam.AlignmentFile(input_bam, "rb") as bam_in:
        for read in bam_in:
            if read.has_tag('SC'):
                sc_tags.add(read.get_tag('SC'))
    return sc_tags

def subset_by_sc_tag(input_bam, output_bam, sc_value):
    """Subset BAM file based on SC tag."""
    with pysam.AlignmentFile(input_bam, "rb") as bam_in:
        header = bam_in.header
        with pysam.AlignmentFile(output_bam, "wb", header=header) as bam_out:
            for read in bam_in:
                if read.has_tag('SC') and read.get_tag('SC') == sc_value:
                    bam_out.write(read)

def main():
    parser = argparse.ArgumentParser(description="Split a BAM file into separate files based on unique SC tags.")
    parser.add_argument('-i', '--input', type=str, required=True, help="Path to the input BAM file.")
    parser.add_argument('-o', '--output_dir', type=str, help="Directory to save the output BAM files. If not specified, files will be saved in the current directory.")
    
    args = parser.parse_args()
    
    input_bam_path = args.input
    output_dir = args.output_dir if args.output_dir else os.getcwd()
    
    if not os.path.isfile(input_bam_path):
        print(f"Input BAM file '{input_bam_path}' does not exist.")
        return

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print(f"Reading unique SC tags from BAM file '{input_bam_path}'")
    sc_tags = get_unique_sc_tags(input_bam_path)

    for sc_tag in sc_tags:
        output_bam_path = os.path.join(output_dir, f"subset_{sc_tag}.bam")
        print(f"Processing SC tag '{sc_tag}'")
        subset_by_sc_tag(input_bam_path, output_bam_path, sc_tag)
    
    print(f"Subset BAM files saved to {output_dir}")

if __name__ == "__main__":
    main()
