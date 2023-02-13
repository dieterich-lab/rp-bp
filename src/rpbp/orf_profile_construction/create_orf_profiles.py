#! /usr/bin/env python3

"""This is the main script for the first part of the
pipeline, it creates a filtered genome profile and
extract the ORF profiles.

Calls:
    create-base-genome-profile
    extract-metagene-profiles (pbiotools.ribo)
    estimate-metagene-profile-bayes-factors (pbiotools.ribo)
    select-periodic-offsets (pbiotools.ribo)
    extract-orf-profiles
"""

import logging
import sys
import argparse
import shlex
import yaml


import pbiotools.utils.pgrm_utils as pgrm_utils
import pbiotools.misc.logging_utils as logging_utils
import pbiotools.misc.shell_utils as shell_utils
import pbiotools.misc.utils as utils

import rpbp.ribo_utils.utils as ribo_utils
import rpbp.ribo_utils.filenames as filenames

from rpbp.defaults import (
    default_num_cpus,
    default_mem,
    star_executable,
    metagene_options,
)

default_models_base = filenames.get_default_models_base()

logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""This script runs all of the processing necessary to
        produce the signals used for ORF translation prediction. In particular, it creates the
        metagene profiles, selected the periodic fragments and generate the ORF profiles.""",
    )

    parser.add_argument("raw_data", help="The raw data file (fastq[.gz])")

    parser.add_argument("config", help="The (yaml) configuration file")

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

    parser.add_argument(
        "--mem", help="The amount of RAM to request", default=default_mem
    )

    parser.add_argument("--tmp", help="The location for temp files", default=None)

    parser.add_argument("--do-not-call", action="store_true")

    parser.add_argument(
        "--overwrite",
        help="""If this flag is present, existing files
        will be overwritten.""",
        action="store_true",
    )

    parser.add_argument(
        "-k",
        "--keep-intermediate-files",
        help="""If this flag is given,
        then all intermediate files will be kept; otherwise, they will be deleted.
        This feature is implemented piecemeal. If the --do-not-call flag is given,
        then nothing will be deleted.""",
        action="store_true",
    )

    logging_utils.add_logging_options(parser)
    pgrm_utils.add_star_options(parser, star_executable)
    pgrm_utils.add_flexbar_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "[create-orf-profiles]: {}".format(" ".join(sys.argv))
    logger.info(msg)

    config = yaml.load(open(args.config), Loader=yaml.FullLoader)

    # check that all of the necessary programs are callable
    programs = [
        "flexbar",
        args.star_executable,
        "samtools",
        "bowtie2",
        "create-base-genome-profile",
        "remove-multimapping-reads",
        "extract-metagene-profiles",
        "estimate-metagene-profile-bayes-factors",
        "select-periodic-offsets",
        "extract-orf-profiles",
    ]
    shell_utils.check_programs_exist(programs)

    required_keys = [
        "riboseq_data",
        "ribosomal_index",
        "genome_base_path",
        "genome_name",
    ]
    utils.check_keys_exist(config, required_keys)

    note = config.get("note", None)
    models_base = config.get("models_base", default_models_base)

    logging_str = logging_utils.get_logging_options_string(args)
    star_str = pgrm_utils.get_star_options_string(args)
    flexbar_str = pgrm_utils.get_flexbar_options_string(args)

    # handle do_not_call so that we do call the preprocessing script,
    # but that it does not run anything
    call = not args.do_not_call
    do_not_call_argument = ""
    if not call:
        do_not_call_argument = "--do-not-call"

    overwrite_argument = ""
    if args.overwrite:
        overwrite_argument = "--overwrite"

    keep_intermediate_str = ""
    if args.keep_intermediate_files:
        keep_intermediate_str = "--keep-intermediate-files"

    tmp_str = ""
    if args.tmp is not None:
        tmp_str = "--tmp {}".format(args.tmp)

    mem_str = "--mem {}".format(shlex.quote(args.mem))

    # check if we want to keep multimappers
    is_unique = not config.get("keep_riboseq_multimappers", False)

    riboseq_raw_data = args.raw_data
    riboseq_bam_filename = filenames.get_riboseq_bam(
        config["riboseq_data"], args.name, is_unique=is_unique, note=note
    )

    cmd = "create-base-genome-profile {} {} {} --num-cpus {} {} {} {} {} {} {} {} {}".format(
        riboseq_raw_data,
        args.config,
        args.name,
        args.num_cpus,
        do_not_call_argument,
        overwrite_argument,
        logging_str,
        star_str,
        tmp_str,
        flexbar_str,
        keep_intermediate_str,
        mem_str,
    )

    # There could be cases where we start somewhere in the middle of creating
    # the base genome profile. So even if the "raw data" is not available,
    # we still want to call the base pipeline.
    # in_files = [riboseq_raw_data]
    in_files = []
    out_files = [riboseq_bam_filename]
    # we always call this, and pass --do-not-call through
    shell_utils.call_if_not_exists(
        cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=True
    )

    # Extract the metagene profiles

    start_upstream_str = utils.get_config_argument(
        config,
        "metagene_start_upstream",
        "start-upstream",
        default=metagene_options["metagene_start_upstream"],
    )
    start_downstream_str = utils.get_config_argument(
        config,
        "metagene_start_downstream",
        "start-downstream",
        default=metagene_options["metagene_start_downstream"],
    )
    end_upstream_str = utils.get_config_argument(
        config,
        "metagene_end_upstream",
        "end-upstream",
        default=metagene_options["metagene_end_upstream"],
    )
    end_downstream_str = utils.get_config_argument(
        config,
        "metagene_end_downstream",
        "end-downstream",
        default=metagene_options["metagene_end_downstream"],
    )

    metagene_profiles = filenames.get_metagene_profiles(
        config["riboseq_data"], args.name, is_unique=is_unique, note=note
    )

    # use the canonical transcripts for extracting the metagene profiles
    transcript_bed = filenames.get_bed(
        config["genome_base_path"],
        config["genome_name"],
        is_annotated=True,
    )

    cmd = "extract-metagene-profiles {} {} {} --num-cpus {} {} {} {} {} {}".format(
        riboseq_bam_filename,
        transcript_bed,
        metagene_profiles,
        args.num_cpus,
        logging_str,
        start_upstream_str,
        start_downstream_str,
        end_upstream_str,
        end_downstream_str,
    )

    in_files = [riboseq_bam_filename, transcript_bed]
    out_files = [metagene_profiles]
    file_checkers = {metagene_profiles: utils.check_gzip_file}
    shell_utils.call_if_not_exists(
        cmd,
        out_files,
        in_files=in_files,
        file_checkers=file_checkers,
        overwrite=args.overwrite,
        call=call,
    )

    # estimate the periodicity for each offset for all read lengths
    metagene_profile_bayes_factors = filenames.get_metagene_profiles_bayes_factors(
        config["riboseq_data"], args.name, is_unique=is_unique, note=note
    )

    periodic_models = filenames.get_models(models_base, "periodic")
    non_periodic_models = filenames.get_models(models_base, "nonperiodic")

    periodic_models_str = " ".join(periodic_models)
    non_periodic_models_str = " ".join(non_periodic_models)

    periodic_models_str = "--periodic-models {}".format(periodic_models_str)
    non_periodic_models_str = "--nonperiodic-models {}".format(non_periodic_models_str)

    periodic_offset_start_str = utils.get_config_argument(
        config,
        "periodic_offset_start",
        default=metagene_options["periodic_offset_start"],
    )
    periodic_offset_end_str = utils.get_config_argument(
        config, "periodic_offset_end", default=metagene_options["periodic_offset_end"]
    )
    metagene_profile_length_str = utils.get_config_argument(
        config,
        "metagene_profile_length",
        default=metagene_options["metagene_profile_length"],
    )
    seed_str = utils.get_config_argument(
        config, "seed", default=metagene_options["seed"]
    )
    chains_str = utils.get_config_argument(
        config, "chains", default=metagene_options["chains"]
    )
    iterations_str = utils.get_config_argument(
        config,
        "metagene_iterations",
        "iterations",
        default=metagene_options["metagene_iterations"],
    )

    cmd = (
        "estimate-metagene-profile-bayes-factors {} {} --num-cpus {} {} {} "
        "{} {} {} {} {} {} {}".format(
            metagene_profiles,
            metagene_profile_bayes_factors,
            args.num_cpus,
            periodic_models_str,
            non_periodic_models_str,
            periodic_offset_start_str,
            periodic_offset_end_str,
            metagene_profile_length_str,
            seed_str,
            chains_str,
            iterations_str,
            logging_str,
        )
    )

    in_files = [metagene_profiles]
    in_files.extend(periodic_models)
    in_files.extend(non_periodic_models)
    out_files = [metagene_profile_bayes_factors]
    file_checkers = {metagene_profile_bayes_factors: utils.check_gzip_file}
    shell_utils.call_if_not_exists(
        cmd,
        out_files,
        in_files=in_files,
        file_checkers=file_checkers,
        overwrite=args.overwrite,
        call=call,
    )

    # select the best read lengths for constructing the signal
    periodic_offsets = filenames.get_periodic_offsets(
        config["riboseq_data"], args.name, is_unique=is_unique, note=note
    )

    cmd = "select-periodic-offsets {} {}".format(
        metagene_profile_bayes_factors, periodic_offsets
    )

    in_files = [metagene_profile_bayes_factors]
    out_files = [periodic_offsets]
    file_checkers = {periodic_offsets: utils.check_gzip_file}
    shell_utils.call_if_not_exists(
        cmd,
        out_files,
        in_files=in_files,
        file_checkers=file_checkers,
        overwrite=args.overwrite,
        call=call,
    )

    # get the lengths and offsets which meet the required criteria from the config file
    lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(
        config,
        args.name,
        args.do_not_call,
        is_unique=is_unique,
        default_params=metagene_options,
    )

    if len(lengths) == 0:
        msg = (
            "No periodic read lengths and offsets were found. Try relaxing "
            "min_metagene_profile_count, min_metagene_bf_mean, max_metagene_bf_var, "
            "and/or min_metagene_bf_likelihood. Quitting."
        )
        logger.critical(msg)
        return

    lengths_str = " ".join(lengths)
    offsets_str = " ".join(offsets)

    seqname_prefix_str = utils.get_config_argument(config, "seqname_prefix")

    # extract the riboseq profiles for each orf
    unique_filename = filenames.get_riboseq_bam(
        config["riboseq_data"], args.name, is_unique=is_unique, note=note
    )

    profiles_filename = filenames.get_riboseq_profiles(
        config["riboseq_data"],
        args.name,
        length=lengths,
        offset=offsets,
        is_unique=is_unique,
        note=note,
    )

    orfs_genomic = filenames.get_orfs(
        config["genome_base_path"], config["genome_name"], note=config.get("orf_note")
    )

    exons_file = filenames.get_exons(
        config["genome_base_path"], config["genome_name"], note=config.get("orf_note")
    )

    cmd = "extract-orf-profiles {} {} {} {} --lengths {} --offsets {} {} {} --num-cpus {} ".format(
        unique_filename,
        orfs_genomic,
        exons_file,
        profiles_filename,
        lengths_str,
        offsets_str,
        logging_str,
        seqname_prefix_str,
        args.num_cpus,
    )

    in_files = [orfs_genomic, exons_file, unique_filename]
    out_files = [profiles_filename]

    # todo: implement a file checker for mtx files
    shell_utils.call_if_not_exists(
        cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call
    )


if __name__ == "__main__":
    main()
