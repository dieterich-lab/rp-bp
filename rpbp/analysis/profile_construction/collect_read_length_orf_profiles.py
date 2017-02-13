#! /usr/bin/env python3

import argparse
import gzip
import string
import scipy.io

import logging
import misc.logging_utils as logging_utils
logger = logging.getLogger(__name__)

default_lengths = []
default_offsets = []

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Collect the individual read length ORF profiles created "
        "by create-read-length-orf-profiles into a single sparse \"tensor\". "
        "N.B. This script is called by create-read-length-orf-profiles, and "
        "it is unlikely that it should ever be called independently.")

    parser.add_argument('mtx_template', help="The template for the individual "
        "output files. The tempalte should use ${length} and ${offset} to "
        "indicate where the respective values should appear in the filename. "
        "N.B. This string probably needs to be given in single quotes so that "
        "the template variables are not interpreted by the shell.")

    parser.add_argument('out', help="The (mtx.gz) output file containing the "
        "ORF profiles and read lengths")


    parser.add_argument('-l', '--lengths', help="If any values are given, "
        "then only reads which have those lengths will be included in the "
        "signal construction.", type=int, default=default_lengths, nargs='*')
    parser.add_argument('-o', '--offsets', help="The 5' end of reads will be "
        "shifted by this amount. There must be one offset value for each "
        "length (given by the --lengths argument.", type=int, 
        default=default_offsets, nargs='*')

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "Creating mtx template"
    logger.info(msg)
    mtx_template = string.Template(args.mtx_template)

    with gzip.open(args.out, 'wb') as target_gz:

        for length, offset in zip(args.lengths, args.offsets):
            mtx = mtx_template.substitute(length=length, offset=offset)
            mtx = scipy.io.mmread(mtx)
            
            msg = "Processing ORF profiles. length: {}. offset: {}.".format(length, offset)
            logger.info(msg)
            
            for row, col, val in zip(mtx.row, mtx.col, mtx.data):
                s = "{} {} {} {}\n".format(length, row, col, val)
                target_gz.write(s.encode())

if __name__ == '__main__':
    main()
