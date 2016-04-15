#! /usr/bin/env python3

import argparse
import os

import misc.utils as utils

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script removes all reads with multiple alignments from a bam file. "
        "It then sorts the reads. Unless specified otherwise, the sorted file is indexed. "
        "The script attempts to properly handle output file names whether align_out ends "
        "in \".bam\" or not.")
    parser.add_argument('align_in', help="The input alignment file")
    parser.add_argument('align_out', help="The output alignment file with multimappers removed")
    parser.add_argument('--do-not-index', help="If this flag is present, then the index WILL "
        "NOT be created.", action='store_true')
    parser.add_argument('--do-not-call', action='store_true')
    args = parser.parse_args()

    programs = ['samtools']
    utils.check_programs_exist(programs)

    call = not args.do_not_call

    # remove the multimappers and sort the remaining reads
    cmd = "samtools view -h {} | grep '^@\|NH:i:1	' | samtools view -bS  - | samtools sort - -o {}".format(args.align_in, align_out)
    utils.check_call(cmd, call=call)

    # index, unless specied not to
    if not args.do_not_index:
        cmd = "samtools index -b {}".format(args.align_out)
        utils.check_call(cmd, call=call)

if __name__ == '__main__':
    main()
