#! /usr/bin/env python3

import argparse
import gzip
import scipy.io
import yaml

import pbiotools.utils.bed_utils as bed_utils

import logging
import pbiotools.misc.logging_utils as logging_utils

import rpbp.ribo_utils.utils as ribo_utils
import rpbp.ribo_utils.filenames as filenames

from rpbp.defaults import metagene_options

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Collect the individual read length ORF profiles (mtx) created "
        "by 'create-read-length-orf-profiles' into a single 'sparse tensor'. "
        "N.B. This script is called by 'create-read-length-orf-profiles', however"
        "we still call each sample independently for condition, lengths and offsets",
    )

    parser.add_argument("config", help="The (yaml) config file")
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

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "Reading config file"
    logger.info(msg)
    config = yaml.load(open(args.config), Loader=yaml.FullLoader)

    # pull out what we need from the config file
    is_unique = not config.get("keep_riboseq_multimappers", False)
    note = config.get("note", None)

    if args.add_ids:
        orf_note = config.get("orf_note", None)
        orfs_file = filenames.get_orfs(
            config["genome_base_path"], config["genome_name"], note=orf_note
        )
        orfs = bed_utils.read_bed(orfs_file)

    names = [args.name]
    if args.is_condition:
        riboseq_replicates = ribo_utils.get_riboseq_replicates(config)
        names = [n for n in riboseq_replicates[args.name]]

    # keep a map from the lengths to the combined profiles
    length_profile_map = {}

    for name in names:

        msg = "Processing sample: {}".format(name)
        logger.info(msg)

        # get the lengths and offsets which meet the required criteria from the config file
        lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(
            config, name, is_unique=is_unique, default_params=metagene_options
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

            mtx = filenames.get_riboseq_profiles(
                config["riboseq_data"],
                name,
                length=[length],
                offset=[offset],
                is_unique=is_unique,
                note=note,
            )

            mtx = scipy.io.mmread(mtx).tocsr()

            prior_mtx = length_profile_map.get(length, None)

            if prior_mtx is None:
                length_profile_map[length] = mtx
            else:
                length_profile_map[length] = prior_mtx + mtx

    if args.add_ids:
        with gzip.open(args.out, "wb") as target_gz:

            for length, mtx in length_profile_map.items():
                mtx = mtx.tocoo()

                msg = "Writing ORF profiles. length: {}.".format(length)
                logger.info(msg)

                for row, col, val in zip(mtx.row, mtx.col, mtx.data):
                    # orf_num are both zero-based, since we are now using coo
                    orf_id = orfs.loc[orfs["orf_num"] == row]["id"].values[0]
                    s = "{} {} {} {} {}\n".format(row, orf_id, col, length, val)
                    target_gz.write(s.encode())
    else:
        with gzip.open(args.out, "wb") as target_gz:

            for length, mtx in length_profile_map.items():
                mtx = mtx.tocoo()

                msg = "Writing ORF profiles. length: {}.".format(length)
                logger.info(msg)

                for row, col, val in zip(mtx.row, mtx.col, mtx.data):
                    s = "{} {} {} {}\n".format(row, col, length, val)
                    target_gz.write(s.encode())


if __name__ == "__main__":
    main()
