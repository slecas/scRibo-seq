#!/usr/bin/env python 

import argparse
import pysam

def trim_and_extract(read1_file, read2_file, output1_file, output2_file):
    print(f"Processing {read1_file} and {read2_file}...")

    with pysam.FastxFile(read1_file) as read1_in, pysam.FastxFile(read2_file) as read2_in:
        with open(output1_file, 'w') as read1_out, open(output2_file, 'w') as read2_out:
            print(f"Trimming first 10 bases from {read1_file} and integrating into {read2_file}...")

            for entry1, entry2 in zip(read1_in, read2_in):
                # Trim first 10 bases from read1
                trimmed_seq1 = entry1.sequence[10:]
                trimmed_qual1 = entry1.quality[10:]
                
                # Extract first 10 bases and qualities
                first10_bases = entry1.sequence[:10]
                first10_quals = entry1.quality[:10]

                # Write modified read1 to output1
                read1_out.write(f"@{entry1.name}\n{trimmed_seq1}\n+\n{trimmed_qual1}\n")

                # Integrate first 10 bases and qualities into read2
                new_seq2 = first10_bases + entry2.sequence
                new_qual2 = first10_quals + entry2.quality

                # Write modified read2 to output2
                read2_out.write(f"@{entry2.name}\n{new_seq2}\n+\n{new_qual2}\n")

    print("Processing completed.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Trim and extract sequences from FASTQ files.')
    parser.add_argument('read1', type=str, help='Path to the first read FASTQ file.')
    parser.add_argument('read2', type=str, help='Path to the second read FASTQ file.')
    parser.add_argument('output1', type=str, help='Path to the output file for the first read.')
    parser.add_argument('output2', type=str, help='Path to the output file for the second read.')

    args = parser.parse_args()

    trim_and_extract(args.read1, args.read2, args.output1, args.output2)