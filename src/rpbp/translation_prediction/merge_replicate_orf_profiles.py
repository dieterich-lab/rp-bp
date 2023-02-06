#! /usr/bin/env python3

import argparse
import logging
import scipy.io

import pbiotools.misc.logging_utils as logging_utils
import pbiotools.misc.math_utils as math_utils

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""This script adds the ORF profiles from a set
        of profiles (presumably, each file corresponds to one replicate from a condition).
        The script keeps the profiles in sparse matrix format, so it is fairly efficient.""",
    )

    parser.add_argument(
        "profiles", help="The (mtx) files containing the ORF profiles", nargs="+"
    )

    parser.add_argument(
        "out", help="The (mtx.gz) output file containing the merged profiles"
    )

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "Reading first ORF profile"
    logger.info(msg)

    merged_profiles = scipy.io.mmread(args.profiles[0]).tocsr()

    msg = "Adding each additional profile"
    logger.info(msg)

    for profile_file in args.profiles[1:]:
        msg = "Reading file: {}".format(profile_file)
        logger.info(msg)

        profiles = scipy.io.mmread(profile_file).tocsr()
        merged_profiles = merged_profiles + profiles

    msg = "Writing merged profiles to disk"
    logger.info(msg)

    math_utils.write_sparse_matrix(args.out, merged_profiles)


if __name__ == "__main__":
    main()
