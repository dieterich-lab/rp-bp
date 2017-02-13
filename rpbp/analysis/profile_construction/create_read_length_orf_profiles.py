#! /usr/bin/env python3

import argparse
import shlex
import string
import sys

import misc.slurm as slurm
import misc.utils as utils

import logging
import misc.logging_utils as logging_utils
logger = logging.getLogger(__name__)

default_lengths = []
default_offsets = []
default_seqname_prefix = ""

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Extract the ORF profiles for each specified read length "
        "and offset independently. One sparse matrix file will be created for "
        "each read length. It then collects the values into a sparse tensor.")

    parser.add_argument('bam', help="The bam file including filtered (unique, "
        "etc.) alignments")
    parser.add_argument('orfs', help="The (bed12) file containing the ORFs")
    parser.add_argument('exons', help="The (bed6+2) file containing the exons")
    
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
       
    parser.add_argument('--seqname-prefix', help="If present, this string "
        "will be prepended to the seqname field of the ORFs.", 
        default=default_seqname_prefix)

    slurm.add_sbatch_options(parser)
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    logging_str = logging_utils.get_logging_options_string(args)
    seqname_str = "--seqname-prefix {}".format(shlex.quote(args.seqname_prefix))
    cpus_str = "--num-cpus {}".format(args.num_cpus)

    mtx_template = string.Template(args.mtx_template)

    msg = "[create-read-length-orf-profiles]: {}".format(' '.join(sys.argv))
    logger.info(msg)
    
    # make sure the number of lengths and offsets match
    if len(args.lengths) != len(args.offsets):
        msg = ("[create-read-length-orf-profiles]: the number of --lengths "
            "and --offsets do not match.")
        raise ValueError(msg)

    # make sure the necessary files exist
    required_files = [args.bam, args.orfs, args.exons]
    msg = "[create-read-length-orf-profiles]: Some input files were missing: "
    utils.check_files_exist(required_files, msg=msg)

    job_ids = []
    for length, offset in zip(args.lengths, args.offsets):
        lengths_str = "--lengths {}".format(length)
        offsets_str = "--offsets {}".format(offset)

        mtx = mtx_template.substitute(length=length, offset=offset)

        cmd = "extract-orf-profiles {} {} {} {} {} {} {} {}".format(
            args.bam,
            args.orfs,
            args.exons,
            mtx,
            lengths_str,
            offsets_str,
            seqname_str,
            cpus_str,
            logging_str
        )
        
        job_id = slurm.check_sbatch(cmd, args=args)

        job_ids.append(job_id)

    # now, collect them into a single file
    offsets_str = ' '.join(str(o) for o in args.offsets)
    lengths_str = ' '.join(str(l) for l in args.lengths)

    offsets_str = "--offsets {}".format(offsets_str)
    lengths_str = "--lengths {}".format(lengths_str)

    cmd = "collect-read-length-orf-profiles {} {} {} {} {}".format(
        shlex.quote(args.mtx_template),
        args.out,
        offsets_str,
        lengths_str,
        logging_str
    )

    slurm.check_sbatch(cmd, args=args, dependencies=job_ids)

if __name__ == '__main__':
    main()
