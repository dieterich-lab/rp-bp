#! /usr/bin/env python3

import logging
import sys
import argparse
import os
import pandas as pd

import yaml

import misc.utils as utils

import rpbp.filenames as filenames

default_num_procs = 2
default_star_executable = "STAR"


def main():
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script runs all of the processing necessary to produce the "
        "signals used for later processing. In particular, it runs the standard "
        "rnaseq and riboseq preprocessing, estimates the abundance of transcripts with "
        "the rnaseq data, and selects the most-expressed isoforms and ORFs. Next, it removes "
        "multimapping and non-periodic-length riboseq reads. Finally, it extracts the riboseq "
        "signal for the most-expressed transcripts.")
    parser.add_argument('raw_data', help="The raw data file (fastq[.gz])")
    parser.add_argument('config', help="The (json) config file")
    parser.add_argument('name', help="The name for the dataset, used in the created files")

    parser.add_argument('-p', '--num-procs', help="The number of processors to use",
        type=int, default=default_num_procs)

    parser.add_argument('--star-executable', help="The name of the STAR executable",
        default=default_star_executable)

    parser.add_argument('--do-not-call', action='store_true')
    parser.add_argument('--overwrite', help="If this flag is present, existing files "
        "will be overwritten.", action='store_true')
            
    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    logging_str = utils.get_logging_options_string(args)

    config = yaml.load(open(args.config))
    call = not args.do_not_call

    # check that all of the necessary programs are callable
    programs =  [   'flexbar',
                    args.star_executable,
                    'samtools',
                    'bowtie2',
                    'create-base-genome-profile',
                    'remove-multimapping-reads',
                    'filter-alignment-lengths',
                    'extract-metagene-profiles',
                    'estimate-metagene-profile-bayes-factors',
                    'select-periodic-offsets',
                ]
    utils.check_programs_exist(programs)

    
    required_keys = [   'riboseq_data',
                        'ribosomal_index',
                        'gtf',
                        'star_index',
                        'orfs',
                        'periodic_models',
                        'nonperiodic_models'
                    ]
    utils.check_keys_exist(config, required_keys)


    # the first step is the standard riboseq preprocessing
    
    # handle do_not_call so that we _do_ call the preprocessing script, but that it does not run anything
    do_not_call_argument = ""
    if not call:
        do_not_call_argument = "--do-not-call"

    overwrite_argument = ""
    if args.overwrite:
        overwrite_argument = "--overwrite"

    star_str = "--star-executable {}".format(args.star_executable)
    riboseq_raw_data = args.raw_data
    riboseq_bam_filename = filenames.get_riboseq_bam(config['riboseq_data'], args.name, is_unique=True)
    cmd = "create-base-genome-profile {} {} {} --num-procs {} {} {} {} {}".format(riboseq_raw_data, 
        args.config, args.name, args.num_procs, do_not_call_argument, overwrite_argument, logging_str, star_str)
    in_files = [riboseq_raw_data]
    out_files = [riboseq_bam_filename]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=True) # we always call this, and pass --do-not-call through

    # create the metagene profiles
    metagene_profiles = filenames.get_metagene_profiles(config['riboseq_data'], args.name, is_unique=True)

    seqids_to_keep_str = utils.get_config_argument(config, 'seqids_to_keep')
    start_upstream_str = utils.get_config_argument(config, 
        'metagene_profile_start_upstream', 'start-upstream')
    start_downstream_str = utils.get_config_argument(config, 
        'metagene_profile_start_downstream', 'start-downstream')
    end_upstream_str = utils.get_config_argument(config, 
        'metagene_profile_end_upstream', 'end-upstream')
    end_downstream_str = utils.get_config_argument(config, 
        'metagene_profile_end_downstream', 'end-downstream')

    cmd = "extract-metagene-profiles {} {} {} --num-procs {} {} {} {} {} {} {}".format(riboseq_bam_filename,
        config['orfs'], metagene_profiles, args.num_procs, logging_str, seqids_to_keep_str,
        start_upstream_str, start_downstream_str, end_upstream_str, end_downstream_str)

    in_files = [riboseq_bam_filename, config['orfs']]
    out_files = [metagene_profiles]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)

    # estimate the periodicity for each offset for all read lengths
    metagene_profile_bayes_factors = filenames.get_metagene_profiles_bayes_factors(config['riboseq_data'],
        args.name, is_unique=True)

    periodic_models_str = utils.get_config_argument(config, 'periodic_models')
    non_periodic_models_str = utils.get_config_argument(config, 'nonperiodic_models')
    periodic_offset_start_str = utils.get_config_argument(config, 'periodic_offset_start')
    periodic_offset_end_str = utils.get_config_argument(config, 'periodic_offset_end')
    metagene_profile_length_str = utils.get_config_argument(config, 'metagene_profile_length')
    seed_str = utils.get_config_argument(config, 'seed')
    chains_str = utils.get_config_argument(config, 'chains')
    iterations_str = utils.get_config_argument(config, 'metagene_profile_iterations', 'iterations')

    cmd = ("estimate-metagene-profile-bayes-factors {} {} --num-procs {} {} {} "
        "{} {} {} {} {} {} {}".format(metagene_profiles, 
        metagene_profile_bayes_factors, args.num_procs, periodic_models_str, non_periodic_models_str,
        periodic_offset_start_str, periodic_offset_end_str, metagene_profile_length_str,
        seed_str, chains_str, iterations_str, logging_str))

    in_files = [metagene_profiles]
    in_files.extend(config['periodic_models'])
    in_files.extend(config['nonperiodic_models'])
    out_files = [metagene_profile_bayes_factors]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)
    
    # select the best read lengths for constructing the signal
    periodic_offsets = filenames.get_periodic_offsets(config['riboseq_data'], args.name, is_unique=True)

    cmd = "select-periodic-offsets {} {}".format(metagene_profile_bayes_factors, periodic_offsets)
    in_files = [metagene_profile_bayes_factors]
    out_files = [periodic_offsets]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)
   
if __name__ == '__main__':
    main()

