#! /usr/bin/env python3

import logging
import sys
import argparse

import yaml

import pbiotools.misc.logging_utils as logging_utils
import pbiotools.misc.shell_utils as shell_utils
import pbiotools.misc.utils as utils

import rpbp.ribo_utils.utils as ribo_utils
import rpbp.ribo_utils.filenames as filenames

from rpbp.defaults import (
    default_num_cpus,
    translation_options,
    metagene_options,
)

logger = logging.getLogger(__name__)

default_models_base = filenames.get_default_models_base()


def get_profile(name, config, args):
    """This helper function constructs the name of the smooth profile file
    from the given parameters.
    """
    # keep multimappers?
    is_unique = not config.get("keep_riboseq_multimappers", False)

    # get the lengths and offsets which meet the required criteria from the config file
    lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(
        config,
        name,
        args.do_not_call,
        is_unique=is_unique,
        default_params=metagene_options,
    )

    note_str = config.get("note", None)

    if len(lengths) == 0:
        msg = (
            "No periodic read lengths and offsets were found. Try relaxing "
            "min_metagene_profile_count, min_metagene_bf_mean, max_metagene_bf_var, "
            "and/or min_metagene_bf_likelihood."
        )
        logger.critical(msg)
        return

    profiles = filenames.get_riboseq_profiles(
        config["riboseq_data"],
        name,
        length=lengths,
        offset=offsets,
        is_unique=is_unique,
        note=note_str,
    )

    return profiles


def main():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=""""This script runs the second part of the pipeline:
        it estimate ORF Bayes factors using the ORF profiles, then make the final prediction set.""",
    )

    parser.add_argument("config", help="The (yaml) config file")

    parser.add_argument(
        "name", help="The name for the dataset, used in the created files"
    )

    parser.add_argument(
        "-p",
        "--num-cpus",
        help="The number of processors to use",
        type=int,
        default=default_num_cpus,
    )

    parser.add_argument("--do-not-call", action="store_true")

    parser.add_argument(
        "--overwrite",
        help="If this flag is present, existing files will be overwritten.",
        action="store_true",
    )

    parser.add_argument(
        "--merge-replicates",
        help="""If this flag is present, then the ORF profiles
        will be merged for all replicates in the condition given by <name>. The filenames, etc.,
        will reflect the condition name, but not the lengths and offsets of the individual replicates.
        N.B. If this flag is is present, the --overwrite flag will automatically be set!""",
        action="store_true",
    )

    parser.add_argument(
        "--write-unfiltered",
        help="""If this flag is given, in addition to the default
        filtered predictions (longest ORF for each stop codon, then
        highest Bayes factor among overlapping ORFs), output all ORF
        predictions""",
        action="store_true",
    )

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "[predict_translated_orfs]: {}".format(" ".join(sys.argv))
    logger.debug(msg)

    logging_str = logging_utils.get_logging_options_string(args)

    config = yaml.load(open(args.config), Loader=yaml.FullLoader)
    call = not args.do_not_call

    # check that all of the necessary programs are callable
    programs = ["estimate-orf-bayes-factors", "select-final-prediction-set"]
    shell_utils.check_programs_exist(programs)

    required_keys = ["riboseq_data", "fasta", "genome_base_path", "genome_name"]
    utils.check_keys_exist(config, required_keys)

    note_str = config.get("note", None)

    # we always need the ORFs
    orfs_genomic = filenames.get_orfs(
        config["genome_base_path"], config["genome_name"], note=config.get("orf_note")
    )

    # smoothing parameters (only for filenames), get actual arguments below
    # default values are not used in the file names
    fraction_name = config.get("smoothing_fraction", None)
    reweighting_iterations_name = config.get("smoothing_reweighting_iterations", None)

    # check if we are running Rp-Bp (default) or Rp-chi
    chi_square_only_str = ""
    chi_square_only = False
    if "chi_square_only" in config:
        chi_square_only_str = "--chi-square-only"
        chi_square_only = True
        fraction_name = None
        reweighting_iterations_name = None
        msg = """ The final prediction set will be made based on the chi square test only!
                  The translation models will not be fit to the data, and the posterior
                  distributions will not be estimated. """
        logger.info(msg)

    # keep multimappers?
    is_unique = not config.get("keep_riboseq_multimappers", False)

    # first, check if we are merging replicates

    # either way, the following variables need to have values for the rest of
    # the pipeline: lengths, offsets, smooth_profiles
    if args.merge_replicates:
        msg = "The --merge-replicates option was given, so --overwrite is being set to True."
        logger.warning(msg)
        args.overwrite = True

        # now, actually merge the replicates
        riboseq_replicates = ribo_utils.get_riboseq_replicates(config)

        # we will not use the lengths and offsets in the filenames
        lengths = None
        offsets = None

        # we will also merge all of unsmoothed profiles
        replicate_profiles = [
            get_profile(name, config, args) for name in riboseq_replicates[args.name]
        ]

        replicate_profiles_str = " ".join(replicate_profiles)

        profiles = filenames.get_riboseq_profiles(
            config["riboseq_data"],
            args.name,
            length=lengths,
            offset=offsets,
            is_unique=is_unique,
            note=note_str,
        )

        cmd = "merge-replicate-orf-profiles {} {} {}".format(
            replicate_profiles_str, profiles, logging_str
        )

        in_files = replicate_profiles
        out_files = [profiles]

        # todo: implement file checker for mtx files
        shell_utils.call_if_not_exists(
            cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call
        )

    else:
        # otherwise, just treat things as normal
        # get the lengths and offsets which meet the required criteria from
        # the config file
        lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(
            config,
            args.name,
            args.do_not_call,
            is_unique=is_unique,
            default_params=metagene_options,
        )

        profiles = get_profile(args.name, config, args)

    # estimate the bayes factors
    bayes_factors = filenames.get_riboseq_bayes_factors(
        config["riboseq_data"],
        args.name,
        length=lengths,
        offset=offsets,
        is_unique=is_unique,
        note=note_str,
        fraction=fraction_name,
        reweighting_iterations=reweighting_iterations_name,
    )

    # the smoothing options
    min_length_str = utils.get_config_argument(
        config,
        "min_orf_length",
        "min-length",
        default=translation_options["orf_min_length_pre"],
    )

    max_length_str = utils.get_config_argument(
        config,
        "max_orf_length",
        "max-length",
        default=translation_options["orf_max_length_pre"],
    )

    min_profile_str = utils.get_config_argument(
        config,
        "min_signal",
        "min-profile",
        default=translation_options["orf_min_profile_count_pre"],
    )

    fraction_str = utils.get_config_argument(
        config,
        "smoothing_fraction",
        "fraction",
        default=translation_options["smoothing_fraction"],
    )

    reweighting_iterations_str = utils.get_config_argument(
        config,
        "smoothing_reweighting_iterations",
        "reweighting-iterations",
        default=translation_options["smoothing_reweighting_iterations"],
    )

    # parse out all of the options from the config file, if they are present
    models_base = config.get("models_base", default_models_base)
    translated_models = filenames.get_models(models_base, "translated")
    untranslated_models = filenames.get_models(models_base, "untranslated")

    translated_models_str = " ".join(translated_models)
    untranslated_models_str = " ".join(untranslated_models)

    translated_models_str = "--translated-models {}".format(translated_models_str)
    untranslated_models_str = "--untranslated-models {}".format(untranslated_models_str)

    seed_str = utils.get_config_argument(
        config, "seed", default=translation_options["seed"]
    )

    chains_str = utils.get_config_argument(
        config, "chains", default=translation_options["chains"]
    )

    iterations_str = utils.get_config_argument(
        config,
        "translation_iterations",
        "iterations",
        default=translation_options["translation_iterations"],
    )

    cmd = (
        "estimate-orf-bayes-factors {} {} {} {} {} {} {} {} {} {} "
        "{} {} {} {} {} --num-cpus {}".format(
            profiles,
            orfs_genomic,
            bayes_factors,
            translated_models_str,
            untranslated_models_str,
            logging_str,
            min_length_str,
            max_length_str,
            min_profile_str,
            fraction_str,
            reweighting_iterations_str,
            seed_str,
            iterations_str,
            chains_str,
            chi_square_only_str,
            args.num_cpus,
        )
    )

    in_files = [profiles, orfs_genomic]
    in_files.extend(translated_models)
    in_files.extend(untranslated_models)
    out_files = [bayes_factors]
    file_checkers = {bayes_factors: utils.check_gzip_file}
    msg = "estimate-bayes-factors in_files: {}".format(in_files)
    logger.debug(msg)
    shell_utils.call_if_not_exists(
        cmd,
        out_files,
        in_files=in_files,
        file_checkers=file_checkers,
        overwrite=args.overwrite,
        call=call,
    )

    filters = [True]
    if args.write_unfiltered:
        filters.append(False)
    for is_filtered in filters:

        filtered_str = ""
        if is_filtered:
            filtered_str = "--select-longest-by-stop --select-best-overlapping"

        # now, select the ORFs (longest for each stop codon) which pass the prediction filters
        predicted_orfs = filenames.get_riboseq_predicted_orfs(
            config["riboseq_data"],
            args.name,
            length=lengths,
            offset=offsets,
            is_unique=is_unique,
            note=note_str,
            fraction=fraction_name,
            reweighting_iterations=reweighting_iterations_name,
            is_filtered=is_filtered,
            is_chisq=chi_square_only,
        )

        predicted_orfs_dna = filenames.get_riboseq_predicted_orfs_dna(
            config["riboseq_data"],
            args.name,
            length=lengths,
            offset=offsets,
            is_unique=is_unique,
            note=note_str,
            fraction=fraction_name,
            reweighting_iterations=reweighting_iterations_name,
            is_filtered=is_filtered,
            is_chisq=chi_square_only,
        )

        predicted_orfs_protein = filenames.get_riboseq_predicted_orfs_protein(
            config["riboseq_data"],
            args.name,
            length=lengths,
            offset=offsets,
            is_unique=is_unique,
            note=note_str,
            fraction=fraction_name,
            reweighting_iterations=reweighting_iterations_name,
            is_filtered=is_filtered,
            is_chisq=chi_square_only,
        )

        min_bf_mean_str = utils.get_config_argument(
            config, "min_bf_mean", default=translation_options["min_bf_mean"]
        )

        max_bf_var_str = utils.get_config_argument(
            config, "max_bf_var", default=translation_options["max_bf_var"]
        )

        min_bf_likelihood_str = utils.get_config_argument(
            config,
            "min_bf_likelihood",
            default=translation_options["min_bf_likelihood"],
        )

        min_profile_str = utils.get_config_argument(
            config,
            "orf_min_profile_count",
            "min-profile",
            default=translation_options["orf_min_profile_count"],
        )

        min_length_str = utils.get_config_argument(
            config,
            "orf_min_length",
            "min-length",
            default=translation_options["orf_min_length"],
        )

        chisq_significance_level_str = utils.get_config_argument(
            config,
            "chisq_alpha",
            "chisq-significance-level",
            default=translation_options["chisq_alpha"],
        )

        cmd = "select-final-prediction-set {} {} {} {} {} {} {} {} {} {} {} {} {} {}".format(
            bayes_factors,
            config["fasta"],
            predicted_orfs,
            predicted_orfs_dna,
            predicted_orfs_protein,
            min_bf_mean_str,
            max_bf_var_str,
            min_bf_likelihood_str,
            min_profile_str,
            min_length_str,
            logging_str,
            chi_square_only_str,
            chisq_significance_level_str,
            filtered_str,
        )

        in_files = [bayes_factors, config["fasta"]]
        out_files = [predicted_orfs, predicted_orfs_dna, predicted_orfs_protein]

        file_checkers = {predicted_orfs: utils.check_gzip_file}

        # todo: implement file checker for fasta files
        shell_utils.call_if_not_exists(
            cmd,
            out_files,
            in_files=in_files,
            file_checkers=file_checkers,
            overwrite=args.overwrite,
            call=call,
        )


if __name__ == "__main__":
    main()
