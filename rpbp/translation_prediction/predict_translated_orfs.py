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
default_tmp = None # utils.abspath('tmp')

default_min_metagene_profile_count = 1000
default_min_metagene_profile_bayes_factor_mean = 5
default_max_metagene_profile_bayes_factor_var = 5

def get_periodic_lengths_and_offsets(config, name, do_not_call):
    # check if we specified to just use a fixed offset and length
    if 'use_fixed_lengths' in config:
        lengths = config['lengths']
        offsets = config['offsets']

        return (lengths, offsets)

    # filter out the lengths which do not satisfy the quality thresholds
    min_metagene_profile_count = config.get(
        "min_metagene_profile_count", default_min_metagene_profile_count)

    min_metagene_profile_bayes_factor_mean = config.get(
        "min_metagene_profile_bayes_factor_mean", default_min_metagene_profile_bayes_factor_mean)

    max_metagene_profile_bayes_factor_var = config.get(
        "max_metagene_profile_bayes_factor_var", default_max_metagene_profile_bayes_factor_var)

    periodic_offsets = filenames.get_periodic_offsets(config['riboseq_data'], name, is_unique=True)
    
    if not os.path.exists(periodic_offsets) and not do_not_call:
        msg = ("The periodic offsets file does not exist. Please ensure the select-periodic-offsets "
            "script completed successfully or specify the \"use_fixed_lengths\", \"lengths\", and "
            "\"offsets\" values in the configuration file. Filename: {}".format(periodic_offsets))
        raise FileNotFoundError(msg)
    else:
        msg = ("The periodic offsets file does not exist. Please ensure the select-periodic-offsets "
            "script completed successfully or specify the \"use_fixed_lengths\", \"lengths\", and "
            "\"offsets\" values in the configuration file.\n\nThe --do-not-call flag was given, so "
            "\"dummy\" default lengths will be used to check the remaining calls.")
        logging.warning(msg)

        offsets = ["12"]
        lengths = ["29"]
        return (lengths, offsets)
        
    
    offsets_df = pd.read_csv(periodic_offsets)
    m_count = offsets_df['highest_peak_peak'] > min_metagene_profile_count
    m_bf_mean = offsets_df['highest_peak_bf_mean'] > min_metagene_profile_bayes_factor_mean
    m_bf_var = offsets_df['highest_peak_bf_var'] < max_metagene_profile_bayes_factor_var

    filtered_periodic_offsets = offsets_df[m_count & m_bf_mean & m_bf_var]

    offsets = filtered_periodic_offsets['highest_peak_offset']
    lengths = filtered_periodic_offsets['length']

    # offsets must be positive
    offsets = [str(-1*int(o)) for o in offsets]
    lengths = [str(int(l)) for l in lengths]

    return (lengths, offsets)


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
                        'translated_models',
                        'untranslated_models'
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
    lengths, offsets = get_periodic_lengths_and_offsets(config, args.name, args.do_not_call)
    lengths_str = ' '.join(lengths)
    offsets_str = ' '.join(offsets)

    seqname_prefix_str = utils.get_config_argument(config, 'seqname_prefix')
    
    # extract the riboseq profiles for each orf
    unique_filename = filenames.get_riboseq_bam(config['riboseq_data'], args.name, is_unique=True)
    profiles_filename = filenames.get_riboseq_profiles(config['riboseq_data'], args.name, 
        length=lengths, offset=offsets, is_unique=True)

    
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
    translated_models_str = utils.get_config_argument(config, 'translated_models')
    untranslated_models_str = utils.get_config_argument(config, 'untranslated_models')
    
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
    in_files.extend(config['translated_models'])
    in_files.extend(config['untranslated_models'])
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

