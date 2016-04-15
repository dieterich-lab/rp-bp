#! /usr/bin/env python3

import argparse
import os

import misc.utils as utils

default_cds_out = ""

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script extracts the transcript annotations and creates the fasta "
        "file. It requires gffread (from cufflinks).")
    parser.add_argument('gtf', help="The gtf file with the annotations (from ensembl)")
    parser.add_argument('fasta', help="The 'dna.toplevel.fa' fasta file")
    parser.add_argument('out', help="The fasta file containing the transcript sequences. "
        "It will contain the annotated exons (which typically includes UTRs and CDSs.")

    parser.add_argument('--cds-out', help="If this filename is given, then another fasta "
        "file which contains only the CDS sequences is created.", default=default_cds_out)

    args = parser.parse_args()

    # check that all of the necessary programs are callable
    programs =  ['gffread']
    utils.check_programs_exist(programs)

    # first, the entire transcripts
    cmd = "gffread -W -w {} -g {} {}".format(args.out, args.fasta, args.gtf)
    utils.check_call(cmd)

    # now, the cds regions, if desired
    if len(args.cds_out) > 0:
        cmd = "gffread -W -x {} -g {} {}".format(args.cds_out, args.fasta, args.gtf)
        utils.check_call(cmd)

if __name__ == '__main__':
    main()
