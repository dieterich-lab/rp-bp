#! /usr/bin/env python3

import argparse
import sys
import yaml

import pbiotools.misc.slurm as slurm
import pbiotools.misc.utils as utils

import logging
import pbiotools.misc.logging_utils as logging_utils

import rpbp.ribo_utils.utils as ribo_utils
import rpbp.ribo_utils.filenames as filenames

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Extract the ORF profiles for each specified read length "
        "and offset independently, creating one sparse matrix file (mtx) for "
        "each read length. These are then collected into a 'sparse tensor'.",
    )

    parser.add_argument("config", help="The yaml config file.")
    parser.add_argument(
        "name",
        help="The name of either one of the 'riboseq_samples'"
        "or 'riboseq_biological_replicates' from the config file.",
    )

    parser.add_argument(
        "out",
        help="The output (txt.gz) file. N.B. The output uses"
        "base-0 indexing, contrary to the unsmoothed ORF profiles, which are written"
        "using the matrix market format (base-1 indexing).",
    )

    parser.add_argument(
        "-c",
        "--is-condition",
        help="If this flag is present, "
        "then 'name' will be taken to be a condition name. The profiles for "
        "all relevant replicates of the condition will be created.",
        action="store_true",
    )

    parser.add_argument(
        "--add-ids",
        help="If this flag is present, "
        "then orf_ids will be added to the final output.",
        action="store_true",
    )

    slurm.add_sbatch_options(parser)
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    logging_str = logging_utils.get_logging_options_string(args)
    cpus_str = "--num-cpus {}".format(args.num_cpus)

    msg = "[create-read-length-orf-profiles]: {}".format(" ".join(sys.argv))
    logger.info(msg)

    msg = "Reading config file"
    logger.info(msg)
    config = yaml.load(open(args.config), Loader=yaml.FullLoader)

    # pull out what we need from the config file
    is_unique = not config.get("keep_riboseq_multimappers", False)
    seqname_str = utils.get_config_argument(config, "seqname_prefix")
    note = config.get("note", None)
    orf_note = config.get("orf_note", None)

    orfs = filenames.get_orfs(
        config["genome_base_path"], config["genome_name"], note=orf_note
    )

    exons = filenames.get_exons(
        config["genome_base_path"], config["genome_name"], note=orf_note
    )

    # make sure the necessary files exist
    required_files = [orfs, exons]
    msg = "[create-read-length-orf-profiles]: Some input files were missing: "
    utils.check_files_exist(required_files, msg=msg)

    # process one sample or all samples from condition
    names = [args.name]
    is_condition_str = ""
    if args.is_condition:
        is_condition_str = "--is-condition"
        riboseq_replicates = ribo_utils.get_riboseq_replicates(config)
        names = [n for n in riboseq_replicates[args.name]]

    job_ids = []
    for name in names:

        msg = "Processing sample: {}".format(name)
        logger.info(msg)

        # now the relevant files
        bam = filenames.get_riboseq_bam(
            config["riboseq_data"], name, is_unique=is_unique, note=note
        )

        # make sure the necessary files exist
        required_files = [bam]
        msg = "[create-read-length-orf-profiles]: Some input files were missing: "
        utils.check_files_exist(required_files, msg=msg)

        # get the lengths and offsets which meet the required criteria from the config file
        lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(
            config, name, is_unique=is_unique
        )

        if len(lengths) == 0:
            msg = (
                "No periodic read lengths and offsets were found. Try relaxing "
                "min_metagene_profile_count, min_metagene_bf_mean, "
                "max_metagene_bf_var, and/or min_metagene_bf_likelihood. Qutting."
            )
            logger.critical(msg)
            return

        for length, offset in zip(lengths, offsets):
            lengths_str = "--lengths {}".format(length)
            offsets_str = "--offsets {}".format(offset)

            mtx = filenames.get_riboseq_profiles(
                config["riboseq_data"],
                name,
                length=[length],
                offset=[offset],
                is_unique=is_unique,
                note=note,
            )

            cmd = "extract-orf-profiles {} {} {} {} {} {} {} {} {}".format(
                bam,
                orfs,
                exons,
                mtx,
                lengths_str,
                offsets_str,
                seqname_str,
                cpus_str,
                logging_str,
            )

            job_id = slurm.check_sbatch(cmd, args=args)

            job_ids.append(job_id)

    add_ids_str = ""
    if args.add_ids:
        add_ids_str = "--add-ids"

    cmd = "collect-read-length-orf-profiles {} {} {} {} {} {}".format(
        args.config, args.name, args.out, is_condition_str, add_ids_str, logging_str
    )

    slurm.check_sbatch(cmd, args=args, dependencies=job_ids)


if __name__ == "__main__":
    main()
