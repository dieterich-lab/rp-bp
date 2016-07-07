#! /usr/bin/env python3

import logging
import sys
import argparse
import os
import pandas as pd

import yaml

import misc.utils as utils

import riboutils.ribo_utils as ribo_utils
import riboutils.ribo_filenames as filenames

default_num_cpus = 1
default_star_executable = "STAR"
default_tmp = None # utils.abspath("tmp")

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

    parser.add_argument('-p', '--num-cpus', help="The number of processors to use",
        type=int, default=default_num_cpus)

    parser.add_argument('--star-executable', help="The name of the STAR executable",
        default=default_star_executable)
    parser.add_argument('--tmp', help="The location for temp files", default=default_tmp)

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
                    'extract-metagene-profiles',
                    'estimate-metagene-profile-bayes-factors',
                    'select-periodic-offsets',
                    'extract-orf-profiles',
                    'smooth-orf-profiles'
                ]
    utils.check_programs_exist(programs)

    
    required_keys = [   'riboseq_data',
                        'ribosomal_index',
                        'gtf',
                        'genome_base_path',
                        'genome_name',
                        'models_base'
                    ]
    utils.check_keys_exist(config, required_keys)

    note = config.get('note', None)

    star_index = filenames.get_star_index(config['genome_base_path'], 
        config['genome_name'], is_merged=False)

    # the first step is the standard riboseq preprocessing
    
    # handle do_not_call so that we _do_ call the preprocessing script, but that it does not run anything
    do_not_call_argument = ""
    if not call:
        do_not_call_argument = "--do-not-call"

    overwrite_argument = ""
    if args.overwrite:
        overwrite_argument = "--overwrite"

    orfs_genomic = filenames.get_orfs(config['genome_base_path'], config['genome_name'], 
        note=config.get('orf_note'))

    star_str = "--star-executable {}".format(args.star_executable)


    tmp_str = ""
    if args.tmp is not None:
        tmp_str = "--tmp {}".format(args.tmp)

    riboseq_raw_data = args.raw_data
    riboseq_bam_filename = filenames.get_riboseq_bam(config['riboseq_data'], args.name, is_unique=True, note=note)
    cmd = "create-base-genome-profile {} {} {} --num-cpus {} {} {} {} {} {}".format(riboseq_raw_data, 
        args.config, args.name, args.num_cpus, do_not_call_argument, overwrite_argument, 
        logging_str, star_str, tmp_str)

    # There could be cases where we start somewhere in the middle of creating
    # the base genome profile. So even if the "raw data" is not available, 
    # we still want to call the base pipeline.
    #in_files = [riboseq_raw_data]
    in_files = []
    out_files = [riboseq_bam_filename]

    # we always call this, and pass --do-not-call through
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=True) 

    # create the metagene profiles
    metagene_profiles = filenames.get_metagene_profiles(config['riboseq_data'], args.name, is_unique=True, note=note)

    seqids_to_keep_str = utils.get_config_argument(config, 'seqids_to_keep')
    start_upstream_str = utils.get_config_argument(config, 
        'metagene_profile_start_upstream', 'start-upstream')
    start_downstream_str = utils.get_config_argument(config, 
        'metagene_profile_start_downstream', 'start-downstream')
    end_upstream_str = utils.get_config_argument(config, 
        'metagene_profile_end_upstream', 'end-upstream')
    end_downstream_str = utils.get_config_argument(config, 
        'metagene_profile_end_downstream', 'end-downstream')

    cmd = "extract-metagene-profiles {} {} {} --num-cpus {} {} {} {} {} {} {}".format(riboseq_bam_filename,
        orfs_genomic, metagene_profiles, args.num_cpus, logging_str, seqids_to_keep_str,
        start_upstream_str, start_downstream_str, end_upstream_str, end_downstream_str)

    in_files = [riboseq_bam_filename, orfs_genomic]
    out_files = [metagene_profiles]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)

    # estimate the periodicity for each offset for all read lengths
    metagene_profile_bayes_factors = filenames.get_metagene_profiles_bayes_factors(config['riboseq_data'],
        args.name, is_unique=True, note=note)

    #periodic_models_str = utils.get_config_argument(config, 'periodic_models')
    #non_periodic_models_str = utils.get_config_argument(config, 'nonperiodic_models')
    periodic_models = filenames.get_models(config['models_base'], 'periodic')
    non_periodic_models = filenames.get_models(config['models_base'], 'nonperiodic')
    
    periodic_models_str = ' '.join(periodic_models)
    non_periodic_models_str = ' '.join(non_periodic_models)

    periodic_models_str = "--periodic-models {}".format(periodic_models_str)
    non_periodic_models_str = "--nonperiodic-models {}".format(non_periodic_models_str)
    
    periodic_offset_start_str = utils.get_config_argument(config, 'periodic_offset_start')
    periodic_offset_end_str = utils.get_config_argument(config, 'periodic_offset_end')
    metagene_profile_length_str = utils.get_config_argument(config, 'metagene_profile_length')
    seed_str = utils.get_config_argument(config, 'seed')
    chains_str = utils.get_config_argument(config, 'chains')
    iterations_str = utils.get_config_argument(config, 'metagene_profile_iterations', 'iterations')

    cmd = ("estimate-metagene-profile-bayes-factors {} {} --num-cpus {} {} {} "
        "{} {} {} {} {} {} {}".format(metagene_profiles, 
        metagene_profile_bayes_factors, args.num_cpus, periodic_models_str, non_periodic_models_str,
        periodic_offset_start_str, periodic_offset_end_str, metagene_profile_length_str,
        seed_str, chains_str, iterations_str, logging_str))

    in_files = [metagene_profiles]
    in_files.extend(periodic_models)
    in_files.extend(non_periodic_models)
    out_files = [metagene_profile_bayes_factors]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)
    
    # select the best read lengths for constructing the signal
    periodic_offsets = filenames.get_periodic_offsets(config['riboseq_data'], 
        args.name, is_unique=True, note=note)

    cmd = "select-periodic-offsets {} {}".format(metagene_profile_bayes_factors, periodic_offsets)
    in_files = [metagene_profile_bayes_factors]
    out_files = [periodic_offsets]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)

    # get the lengths and offsets which meet the required criteria from the config file
    lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(config, args.name, args.do_not_call)

    if len(lengths) == 0:
        msg = ("No periodic read lengths and offsets were found. Try relaxing "
            "min_metagene_profile_count, min_metagene_bf_mean, max_metagene_bf_var, "
            "and/or min_metagene_bf_likelihood. Qutting.")
        logging.critical(msg)
        return

    lengths_str = ' '.join(lengths)
    offsets_str = ' '.join(offsets)

    seqname_prefix_str = utils.get_config_argument(config, 'seqname_prefix')
    
    # extract the riboseq profiles for each orf
    unique_filename = filenames.get_riboseq_bam(config['riboseq_data'], args.name, 
        is_unique=True, note=note)

    profiles_filename = filenames.get_riboseq_profiles(config['riboseq_data'], args.name, 
        length=lengths, offset=offsets, is_unique=True, note=note)

    
    orfs_genomic = filenames.get_orfs(config['genome_base_path'], config['genome_name'], 
        note=config.get('orf_note'))

    cmd = ("extract-orf-profiles {} {} {} --lengths {} --offsets {} {} {} --num-cpus {} "
        "{}".format(unique_filename, orfs_genomic, profiles_filename, lengths_str, 
            offsets_str, logging_str, seqname_prefix_str, args.num_cpus, tmp_str))
    in_files = [orfs_genomic, unique_filename]
    out_files = [profiles_filename]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)

    # now, smooth the ORF signal
    min_length_str = utils.get_config_argument(config, 'min_orf_length', 'min-length')
    max_length_str = utils.get_config_argument(config, 'max_orf_length', 'max-length')
    min_signal_str = utils.get_config_argument(config, 'min_signal')

    fraction_str = utils.get_config_argument(config, 'smoothing_fraction', 'fraction')
    reweighting_iterations_str = utils.get_config_argument(config, 
        'smoothing_reweighting_iterations', 'reweighting-iterations')

    # we also need the values
    fraction = config.get('smoothing_fraction', None)
    reweighting_iterations = config.get('smoothing_reweighting_iterations', None)

    smooth_profiles = filenames.get_riboseq_profiles(config['riboseq_data'], args.name, 
        length=lengths, offset=offsets, is_unique=True, note=note, is_smooth=True, 
        fraction=fraction, reweighting_iterations=reweighting_iterations)

    cmd = "smooth-orf-profiles {} {} {} {} {} {} {} {} {} --num-cpus {}".format(orfs_genomic, 
        profiles_filename, smooth_profiles, fraction_str, reweighting_iterations_str, 
        logging_str, min_signal_str, min_length_str, max_length_str, args.num_cpus)
    in_files = [orfs_genomic, profiles_filename]
    out_files = [smooth_profiles]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)
   
if __name__ == '__main__':
    main()

