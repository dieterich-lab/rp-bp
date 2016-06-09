#! /usr/bin/env python3

import logging
import sys
import argparse
import os
import pandas as pd

import yaml

import misc.utils as utils

import rpbp.rpbp_utils
import rpbp.filenames as filenames

default_num_procs = 2
default_tmp = None # utils.abspath('tmp')

default_min_metagene_profile_count = 1000
default_min_metagene_profile_bayes_factor_mean = 5
default_max_metagene_profile_bayes_factor_var = 5

def main():
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script runs all of the processing necessary to produce the "
        "signals used for later processing. In particular, it runs the standard "
        "rnaseq and riboseq preprocessing, estimates the abundance of transcripts with "
        "the rnaseq data, and selects the most-expressed isoforms and ORFs. Next, it removes "
        "multimapping and non-periodic-length riboseq reads. Finally, it extracts the riboseq "
        "signal for the most-expressed transcripts.")
    parser.add_argument('config', help="The (yaml) config file")
    parser.add_argument('name', help="The name for the dataset, used in the created files")

    parser.add_argument('-p', '--num-procs', help="The number of processors to use",
        type=int, default=default_num_procs)
    parser.add_argument('--tmp', help="The temp directory for pybedtools", default=default_tmp)
    
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
    programs =  [   
                    'extract-orf-profiles',
                    'estimate-orf-bayes-factors',
                    'select-final-prediction-set',
                    'bamToBed',
                    'fastaFromBed'
                ]
    utils.check_programs_exist(programs)

    
    required_keys = [   'riboseq_data',
                        'fasta',
                        'genome_base_path',
                        'genome_name',
                        'models_base'
                        #'translated_models',
                        #'untranslated_models'
                    ]
    utils.check_keys_exist(config, required_keys)

    note_str = config.get('note', None)


    # the first step is the standard riboseq preprocessing
    
    # handle do_not_call so that we _do_ call the preprocessing script, but that it does not run anything
    do_not_call_argument = ""
    if not call:
        do_not_call_argument = "--do-not-call"

    overwrite_argument = ""
    if args.overwrite:
        overwrite_argument = "--overwrite"

    # get the lengths and offsets which meet the required criteria from the config file
    lengths, offsets = rpbp.rpbp_utils.get_periodic_lengths_and_offsets(config, args.name, args.do_not_call)

    if len(lengths) == 0:
        msg = ("No periodic read lengths and offsets were found. Try relaxing "
            "min_metagene_profile_count, min_metagene_profile_bayes_factor_mean, "
            "max_metagene_profile_bayes_factor_var. Qutting.")
        logging.critical(msg)
        return

    lengths_str = ' '.join(lengths)
    offsets_str = ' '.join(offsets)

    seqname_prefix_str = utils.get_config_argument(config, 'seqname_prefix')
    
    # extract the riboseq profiles for each orf
    unique_filename = filenames.get_riboseq_bam(config['riboseq_data'], args.name, is_unique=True, note=note_str)
    profiles_filename = filenames.get_riboseq_profiles(config['riboseq_data'], args.name, 
        length=lengths, offset=offsets, is_unique=True, note=note_str)

    
    orfs_genomic = filenames.get_orfs(config['genome_base_path'], config['genome_name'], 
        note=config.get('orf_note'))

    cmd = ("extract-orf-profiles {} {} {} --lengths {} --offsets {} {} {} --num-procs {} "
            "--tmp {}".format(unique_filename, orfs_genomic, profiles_filename, lengths_str, 
            offsets_str, logging_str, seqname_prefix_str, args.num_procs, args.tmp))
    in_files = [orfs_genomic, unique_filename]
    out_files = [profiles_filename]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)

    # estimate the bayes factors
    bayes_factors = filenames.get_riboseq_bayes_factors(config['riboseq_data'], args.name, 
        length=lengths, offset=offsets, is_unique=True, note=note_str)
    
    # parse out all of the options from the config file, if they are present
    #translated_models_str = utils.get_config_argument(config, 'translated_models')
    #untranslated_models_str = utils.get_config_argument(config, 'untranslated_models')
    translated_models = filenames.get_models(config['models_base'], 'translated')
    untranslated_models = filenames.get_models(config['models_base'], 'untranslated')

    translated_models_str = ' '.join(translated_models)
    untranslated_models_str = ' '.join(untranslated_models)

    translated_models_str = "--translated-models {}".format(translated_models_str)
    untranslated_models_str = "--untranslated-models {}".format(untranslated_models_str)
    
    orf_types_str = utils.get_config_argument(config, 'orf_types')
    min_length_str = utils.get_config_argument(config, 'min_orf_length', 'min-length')
    max_length_str = utils.get_config_argument(config, 'max_orf_length', 'max-length')
    min_signal_str = utils.get_config_argument(config, 'min_signal')

    seed_str = utils.get_config_argument(config, 'seed')
    chains_str = utils.get_config_argument(config, 'chains', 'chains')
    iterations_str = utils.get_config_argument(config, 'translation_iterations', 'iterations')

    chi_square_only_str = ""
    chi_square_only = False
    if 'chi_square_only' in config:
        chi_square_only = True
        chi_square_only_str = "--chi-square-only"

    cmd = "estimate-orf-bayes-factors {} {} {} {} {} {} {} {} {} {} {} {} {} {} --num-procs {}".format(
        profiles_filename, orfs_genomic, bayes_factors, translated_models_str, untranslated_models_str, 
        logging_str, orf_types_str, min_signal_str, min_length_str, max_length_str, seed_str, 
        iterations_str, chains_str, chi_square_only_str, args.num_procs)
    
    in_files = [profiles_filename, orfs_genomic]
    in_files.extend(translated_models)
    in_files.extend(untranslated_models)
    out_files = [bayes_factors]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)

    # now, select the ORFs (longest for each stop codon) which pass the prediction filters
    predicted_orfs = filenames.get_riboseq_predicted_orfs(config['riboseq_data'], 
        args.name, length=lengths, offset=offsets, is_unique=True, note=note_str)
    predicted_orfs_dna = filenames.get_riboseq_predicted_orfs_dna(config['riboseq_data'], 
        args.name, length=lengths, offset=offsets, is_unique=True, note=note_str)
    predicted_orfs_protein = filenames.get_riboseq_predicted_orfs_protein(
        config['riboseq_data'], args.name, length=lengths, offset=offsets, 
        is_unique=True, note=note_str)

    min_bf_mean_str = utils.get_config_argument(config, 'min_bf_mean')
    max_bf_var_str = utils.get_config_argument(config, 'max_bf_var')
    min_profile_str = utils.get_config_argument(config, 'min_signal', 'minimum-profile-sum')

    cmd = "select-final-prediction-set {} {} {} {} {} {} {} {} {}".format(bayes_factors, 
        config['fasta'], predicted_orfs, predicted_orfs_dna, 
        predicted_orfs_protein, min_bf_mean_str, max_bf_var_str, min_profile_str, logging_str)
    in_files = [bayes_factors, config['fasta']]
    out_files = [predicted_orfs, predicted_orfs_dna, predicted_orfs_protein]
    
    # do not make the call if we want chisq only
    predict_call = call and not chi_square_only
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=predict_call)

    
    # finally, repeat the selection step for chi-square
    predicted_orfs = filenames.get_riboseq_predicted_orfs(config['riboseq_data'], 
        args.name, length=lengths, offset=offsets, is_unique=True, note=note_str, is_chisq=True)
    predicted_orfs_dna = filenames.get_riboseq_predicted_orfs_dna(config['riboseq_data'], 
        args.name, length=lengths, offset=offsets, is_unique=True, note=note_str, is_chisq=True)
    predicted_orfs_protein = filenames.get_riboseq_predicted_orfs_protein(
        config['riboseq_data'], args.name, length=lengths, offset=offsets, 
        is_unique=True, note=note_str, is_chisq=True)

    chisq_significance_level_str = utils.get_config_argument(config, 'chisq_significance_level')
    min_profile_str = utils.get_config_argument(config, 'min_signal', 'minimum-profile-sum')

    cmd = "select-final-prediction-set {} {} {} {} {} {} {} {} --use-chi-square".format(bayes_factors, 
        config['fasta'], predicted_orfs, predicted_orfs_dna, 
        predicted_orfs_protein, chisq_significance_level_str, min_profile_str, logging_str)
    in_files = [bayes_factors, config['fasta']]
    out_files = [predicted_orfs, predicted_orfs_dna, predicted_orfs_protein]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)

if __name__ == '__main__':
    main()

