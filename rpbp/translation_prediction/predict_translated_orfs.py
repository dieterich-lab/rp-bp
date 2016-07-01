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

logger = logging.getLogger(__name__)

default_num_cpus = 2
default_tmp = None # utils.abspath('tmp')

def get_profile(name, is_smooth, config, args):
    """ This helper function constructs the name of the smooth profile file
        from the given parameters.
    """
    # get the lengths and offsets which meet the required criteria from the config file
    lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(config, 
        name, args.do_not_call)

    note_str = config.get('note', None)

    # only use the smoothing options if we are actually using smoothing
    if is_smooth:
        fraction = config.get('smoothing_fraction', None)
        reweighting_iterations = config.get('smoothing_reweighting_iterations', None)
    else:
        fraction = None
        reweighting_iterations = None

    if len(lengths) == 0:
        msg = ("No periodic read lengths and offsets were found. Try relaxing "
            "min_metagene_profile_count, min_metagene_bf_mean, max_metagene_bf_var, "
            "and/or min_metagene_bf_likelihood. Qutting.")
        logger.critical(msg)
        return

    smooth_profiles = filenames.get_riboseq_profiles(config['riboseq_data'], name, 
        length=lengths, offset=offsets, is_unique=True, note=note_str, is_smooth=is_smooth, 
        fraction=fraction, reweighting_iterations=reweighting_iterations)

    return smooth_profiles


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

    parser.add_argument('-p', '--num-cpus', help="The number of processors to use",
        type=int, default=default_num_cpus)
    parser.add_argument('--tmp', help="The temp directory for pybedtools", default=default_tmp)
    
    parser.add_argument('--do-not-call', action='store_true')
    parser.add_argument('--overwrite', help="If this flag is present, existing files "
        "will be overwritten.", action='store_true')

    parser.add_argument('--merge-replicates', help="If this flag is present, then the ORF "
        "profiles will be merged for all replicates in the condition given by <name>. The "
        "filenames, etc., will reflect the condition name, but not the lengths and offsets "
        "of the individual replicates.\n\nN.B. If this flag is is present, the --overwrite "
        "flag will automatically be set!", action='store_true')
        
    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    logging_str = utils.get_logging_options_string(args)

    config = yaml.load(open(args.config))
    call = not args.do_not_call

    # check that all of the necessary programs are callable
    programs =  [   
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
                    ]
    utils.check_keys_exist(config, required_keys)

    note_str = config.get('note', None)

    # we always need the ORFs
    orfs_genomic = filenames.get_orfs(config['genome_base_path'], config['genome_name'], 
        note=config.get('orf_note'))

    # and the smoothing parameters
    fraction = config.get('smoothing_fraction', None)
    reweighting_iterations = config.get('smoothing_reweighting_iterations', None)

    # first, check if we are merging replicates

    # either way, the following variables need to have values for the rest of
    # the pipeline: lengths, offsets, smooth_profiles
    if args.merge_replicates:
        msg = ("The --merge-replicates option was given, so --overwrite is "
            "being set to True.")
        logger.warning(msg)
        args.overwrite = True

        # now, actually merge the replicates
        riboseq_replicates = ribo_utils.get_riboseq_replicates(config)

        # we need all of the replicate smooth profiles
        is_smooth = True
        smooth_replicate_profiles = [
            get_profile(name, is_smooth, config, args) for name in riboseq_replicates[args.name]
        ]

        smooth_replicate_profiles_str = ' '.join(smooth_replicate_profiles)


        # we will not use the lengths and offsets in the filenames
        lengths = None
        offsets = None

        # we will call the merged profiles 'smooth_profiles' to match the rest of the pipeline
        smooth_profiles = filenames.get_riboseq_profiles(config['riboseq_data'], args.name, 
            length=lengths, offset=offsets, is_unique=True, note=note_str, is_smooth=is_smooth, 
            fraction=fraction, reweighting_iterations=reweighting_iterations)

        cmd = "merge-replicate-orf-profiles {} {} {}".format(smooth_replicate_profiles_str,
            smooth_profiles, logging_str)
        in_files = smooth_replicate_profiles
        out_files = [smooth_profiles]
        utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)

        # we will also merge all of unsmoothed profiles
        
        is_smooth = False
        unsmoothed_replicate_profiles = [
            get_profile(name, is_smooth, config, args) for name in riboseq_replicates[args.name]
        ]

        unsmoothed_replicate_profiles_str = ' '.join(unsmoothed_replicate_profiles)

        # we will call the merged profiles 'smooth_profiles' to match the rest of the pipeline
        unsmoothed_profiles = filenames.get_riboseq_profiles(config['riboseq_data'], args.name, 
            length=lengths, offset=offsets, is_unique=True, note=note_str, is_smooth=is_smooth, 
            fraction=fraction, reweighting_iterations=reweighting_iterations)

        cmd = "merge-replicate-orf-profiles {} {} {}".format(unsmoothed_replicate_profiles_str,
            unsmoothed_profiles, logging_str)
        in_files = unsmoothed_replicate_profiles
        out_files = [unsmoothed_profiles]
        utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)


    else:
        # otherwise, just treat things as normal
        # get the lengths and offsets which meet the required criteria from the config file
        lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(config, 
            args.name, args.do_not_call)
        is_smooth = True
        smooth_profiles = get_profile(args.name, is_smooth, config, args)
        is_smooth = False
        unsmoothed_profiles = get_profile(args.name, is_smooth, config, args)
        
    # estimate the bayes factors
    bayes_factors = filenames.get_riboseq_bayes_factors(config['riboseq_data'], args.name, 
        length=lengths, offset=offsets, is_unique=True, note=note_str, is_smooth=True, 
        fraction=fraction, reweighting_iterations=reweighting_iterations)
    
    # parse out all of the options from the config file, if they are present
    translated_models = filenames.get_models(config['models_base'], 'translated')
    untranslated_models = filenames.get_models(config['models_base'], 'untranslated')

    translated_models_str = ' '.join(translated_models)
    untranslated_models_str = ' '.join(untranslated_models)

    translated_models_str = "--translated-models {}".format(translated_models_str)
    untranslated_models_str = "--untranslated-models {}".format(untranslated_models_str)
    
    orf_types_str = utils.get_config_argument(config, 'orf_types')
    
    seed_str = utils.get_config_argument(config, 'seed')
    chains_str = utils.get_config_argument(config, 'chains', 'chains')
    iterations_str = utils.get_config_argument(config, 'translation_iterations', 'iterations')

    chi_square_only_str = ""
    chi_square_only = False
    if 'chi_square_only' in config:
        chi_square_only = True
        chi_square_only_str = "--chi-square-only"

    cmd = "estimate-orf-bayes-factors {} {} {} {} {} {} {} {} {} {} {} --num-cpus {}".format(
        smooth_profiles, orfs_genomic, bayes_factors, translated_models_str, untranslated_models_str, 
        logging_str, orf_types_str, seed_str, 
        iterations_str, chains_str, chi_square_only_str, args.num_cpus)
    
    in_files = [smooth_profiles, orfs_genomic]
    in_files.extend(translated_models)
    in_files.extend(untranslated_models)
    out_files = [bayes_factors]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)

    # now, call this again, but use the unsmoothed data, and only use the chisquare test

    # also, do not include the smoothing variables
    rpchi_pvalues = filenames.get_riboseq_bayes_factors(config['riboseq_data'], args.name, 
        length=lengths, offset=offsets, is_unique=True, note=note_str, is_smooth=False)
    chi_square_only_str = "--chi-square-only"

    cmd = "estimate-orf-bayes-factors {} {} {} {} {} {} {} {} {} {} {} --num-cpus {}".format(
        unsmoothed_profiles, orfs_genomic, rpchi_pvalues, translated_models_str, untranslated_models_str, 
        logging_str, orf_types_str, seed_str, 
        iterations_str, chains_str, chi_square_only_str, args.num_cpus)

    in_files = [unsmoothed_profiles, orfs_genomic]
    in_files.extend(translated_models)
    in_files.extend(untranslated_models)
    out_files = [rpchi_pvalues]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)



    # now, select the ORFs (longest for each stop codon) which pass the prediction filters
    predicted_orfs = filenames.get_riboseq_predicted_orfs(config['riboseq_data'], 
        args.name, length=lengths, offset=offsets, is_unique=True, note=note_str, is_smooth=True, 
        fraction=fraction, reweighting_iterations=reweighting_iterations)
    predicted_orfs_dna = filenames.get_riboseq_predicted_orfs_dna(config['riboseq_data'], 
        args.name, length=lengths, offset=offsets, is_unique=True, note=note_str, is_smooth=True, 
        fraction=fraction, reweighting_iterations=reweighting_iterations)
    predicted_orfs_protein = filenames.get_riboseq_predicted_orfs_protein(
        config['riboseq_data'], args.name, length=lengths, offset=offsets, 
        is_unique=True, note=note_str, is_smooth=True, 
        fraction=fraction, reweighting_iterations=reweighting_iterations)

    min_bf_mean_str = utils.get_config_argument(config, 'min_bf_mean')
    max_bf_var_str = utils.get_config_argument(config, 'max_bf_var')
    min_bf_likelihood_str = utils.get_config_argument(config, 'min_bf_likelihood')

    cmd = "select-final-prediction-set {} {} {} {} {} {} {} {} {}".format(bayes_factors, 
        config['fasta'], predicted_orfs, predicted_orfs_dna, 
        predicted_orfs_protein, min_bf_mean_str, max_bf_var_str, 
        min_bf_likelihood_str, logging_str)
    in_files = [bayes_factors, config['fasta']]
    out_files = [predicted_orfs, predicted_orfs_dna, predicted_orfs_protein]
    
    # do not make the call if we want chisq only
    predict_call = call and not chi_square_only
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=predict_call)
    
    # finally, repeat the selection step for chi-square

    # do not include the smoothing options
    predicted_orfs = filenames.get_riboseq_predicted_orfs(config['riboseq_data'], 
        args.name, length=lengths, offset=offsets, is_unique=True, note=note_str, 
        is_chisq=True, is_smooth=False)
    predicted_orfs_dna = filenames.get_riboseq_predicted_orfs_dna(config['riboseq_data'], 
        args.name, length=lengths, offset=offsets, is_unique=True, note=note_str, 
        is_chisq=True, is_smooth=False)
    predicted_orfs_protein = filenames.get_riboseq_predicted_orfs_protein(
        config['riboseq_data'], args.name, length=lengths, offset=offsets, 
        is_unique=True, note=note_str, is_chisq=True, is_smooth=False)

    chisq_significance_level_str = utils.get_config_argument(config, 'chisq_significance_level')
    min_profile_str = utils.get_config_argument(config, 'min_signal', 'minimum-profile-sum')

    cmd = "select-final-prediction-set {} {} {} {} {} {} {} --use-chi-square".format(rpchi_pvalues, 
        config['fasta'], predicted_orfs, predicted_orfs_dna, 
        predicted_orfs_protein, chisq_significance_level_str, logging_str)
    in_files = [rpchi_pvalues, config['fasta']]
    out_files = [predicted_orfs, predicted_orfs_dna, predicted_orfs_protein]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)

if __name__ == '__main__':
    main()

