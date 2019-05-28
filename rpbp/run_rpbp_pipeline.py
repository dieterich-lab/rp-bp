#! /usr/bin/env python3

import argparse
import logging
import os
import shlex
import sys

import yaml

import pbio.utils.star_utils as star_utils
import pbio.misc.logging_utils as logging_utils
import pbio.misc.shell_utils as shell_utils
import pbio.misc.slurm as slurm
import pbio.misc.utils as utils

logger = logging.getLogger(__name__)

default_num_procs = 1
default_tmp = None  # utils.abspath('tmp')


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""This script runs the Rp-Bp and Rp-chi 
                                     pipelines on a given sample. It requires a YAML config 
                                     file that includes a number of keys. Please see the 
                                     documentation for a complete description.""")

    parser.add_argument('raw_data', help="The raw data file (fastq[.gz])")
    parser.add_argument('config', help="The (yaml) config file")
    parser.add_argument('name', help="The name for the dataset, used in the created files")

    parser.add_argument('--tmp', help="The temp directory", default=default_tmp)

    parser.add_argument('--flexbar-options', help="""Optional argument: a space-delimited 
        list of options to pass to flexbar. Each option must be quoted separately as in 
        "--flexbarOption value", using hard, then soft quotes, where "--flexbarOption" 
        is the long parameter name from flexbar and "value" is the value given to this parameter. 
        If specified, flexbar options will override default settings.""", nargs='*', type=str)
    
    parser.add_argument('--overwrite', help="If this flag is present, existing files "
                                            "will be overwritten.", action='store_true')

    parser.add_argument('--profiles-only', help="If this flag is present, then only "
                                                "the ORF profiles will be created", action='store_true')
         
    parser.add_argument('-k', '--keep-intermediate-files', help="""If this flag is given,
        then all intermediate files will be kept; otherwise, they will be
        deleted. This feature is implemented piecemeal. If the --do-not-call flag
        is given, then nothing will be deleted.""", action='store_true')
           
    star_utils.add_star_options(parser)
    slurm.add_sbatch_options(parser)
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    logging_str = logging_utils.get_logging_options_string(args)
    star_str = star_utils.get_star_options_string(args)

    config = yaml.load(open(args.config))
    call = not args.do_not_call

    # check that all of the necessary programs are callable
    programs = [
        'flexbar',
        args.star_executable,
        'samtools',
        'bowtie2',
        'create-base-genome-profile',
        'remove-multimapping-reads',
        'extract-metagene-profiles',
        'estimate-metagene-profile-bayes-factors',
        'select-periodic-offsets',
        'extract-orf-profiles',
        'estimate-orf-bayes-factors',
        'select-final-prediction-set',
        'create-orf-profiles',
        'predict-translated-orfs'
    ]
    shell_utils.check_programs_exist(programs)

    required_keys = [
        'riboseq_data',
        'ribosomal_index',
        'star_index',
        'genome_base_path',
        'genome_name',
        'fasta',
        'gtf'
    ]
    utils.check_keys_exist(config, required_keys)

    # now, check if we want to use slurm
    msg = "use_slurm: {}".format(args.use_slurm)
    logger.debug(msg)

    # if using slurm, submit the script, but we cannot use sys.argv directly
    # as the shell strips the quotes around the arguments
    if args.use_slurm:
        cmd = "{}".format(' '.join("'" + s + "'" if '"' in s else s for s in sys.argv))
        slurm.check_sbatch(cmd, args=args)
        return

    note_str = config.get('note', None)

    # the first step is the standard riboseq preprocessing
    
    # handle do_not_call so that we _do_ call the preprocessing script, 
    # but that it does not run anything
    do_not_call_str = ""
    if not call:
        do_not_call_str = "--do-not-call"

    overwrite_str = ""
    if args.overwrite:
        overwrite_str = "--overwrite"

    # for a sample, we first create its filtered genome profile
    
    keep_intermediate_str = ""
    if args.keep_intermediate_files:
        keep_intermediate_str = "--keep-intermediate-files"

    tmp_str = ""
    if args.tmp is not None:
        tmp_str = "--tmp {}".format(args.tmp)

    flexbar_option_str = ""
    if args.flexbar_options is not None:
        flexbar_option_str = "--flexbar-options {}".format(' '.join("'" + flx_op + "'"
                                                                    for flx_op in args.flexbar_options))

    mem_str = "--mem {}".format(shlex.quote(args.mem))

    cmd = ("create-orf-profiles {} {} {} --num-cpus {} {} {} {} {} {} {} {} {}".format(
        args.raw_data, args.config, args.name, args.num_cpus, mem_str, do_not_call_str,
        overwrite_str, logging_str, star_str, tmp_str, flexbar_option_str, keep_intermediate_str))

    shell_utils.check_call(cmd)

    # check if we only want to create the profiles
    if args.profiles_only:
        return

    # then we predict the ORFs
    cmd = ("predict-translated-orfs {} {} --num-cpus {} {} {} {}".format(args.config,
                                                                         args.name,
                                                                         args.num_cpus,
                                                                         do_not_call_str,
                                                                         overwrite_str,
                                                                         logging_str))
    shell_utils.check_call(cmd)


if __name__ == '__main__':
    main()
