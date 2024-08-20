#!/usr/bin/env python3

import pysam
import argparse

def create_cr_to_cb_lookup(bam_file):
    cr_to_cb = {}
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            cr_tag = None
            cb_tag = None
            for tag, value in read.tags:
                if tag == 'CR':
                    cr_tag = value
                elif tag == 'CB':
                    cb_tag = value
            if cr_tag and cb_tag:
                cr_to_cb[cr_tag] = cb_tag
    return cr_to_cb

def add_cb_tags(input_bam, output_bam, cr_to_cb):
    with pysam.AlignmentFile(input_bam, "rb") as bam_in:
        with pysam.AlignmentFile(output_bam, "wb", header=bam_in.header) as bam_out:
            for read in bam_in:
                if read.has_tag('CR'):
                    cr_tag = read.get_tag('CR')
                    if cr_tag in cr_to_cb:
                        cb_tag = cr_to_cb[cr_tag]
                        if cb_tag != '-':  
                            read.set_tag('CB', cb_tag, value_type='Z')
                            bam_out.write(read)
                        else:
                            continue
                else:
                    bam_out.write(read)

def main():
    parser = argparse.ArgumentParser(description="Process a BAM file to add CB tags based on CR tags using a lookup from another BAM file.")
    parser.add_argument('lookup_bam', type=str, help="Path to the BAM file with CR and CB tags to create the lookup.")
    parser.add_argument('input_bam', type=str, help="Path to the input BAM file without CB tags.")
    parser.add_argument('output_bam', type=str, help="Path to the output BAM file with added CB tags.")
    
    args = parser.parse_args()

    cr_to_cb = create_cr_to_cb_lookup(args.lookup_bam)

    add_cb_tags(args.input_bam, args.output_bam, cr_to_cb)

if __name__ == "__main__":
    main()  