#! /usr/bin/env python3

"""Extract metagene profiles.

Functions:
    get_interval_df
    get_length_strand_profiles
    get_metagene_profile_df
"""

import argparse
import collections
import numpy as np
import pandas as pd
import sys

import pbiotools.utils.bam_utils as bam_utils
import pbiotools.utils.bed_utils as bed_utils
import pbiotools.misc.pandas_utils as pandas_utils

import logging
import pbiotools.misc.logging_utils as logging_utils

logger = logging.getLogger(__name__)

# When called from the Rp-Bp pipeline (create-orf-profiles), default
# options (or else specified via the configuration file) are always
# passed as arguments.

default_start_upstream = 300
default_start_downstream = 300
default_end_upstream = 300
default_end_downstream = 300


def get_interval_df(start, end, seqname, strand):

    interval_df = pd.DataFrame()
    interval_df["start"] = start
    interval_df["end"] = end
    interval_df["seqname"] = seqname
    interval_df["strand"] = strand
    interval_df["id"] = "."
    interval_df["score"] = 0

    return interval_df


def get_length_strand_profiles(matches, profile_length):

    init = lambda: np.zeros(profile_length, int)
    length_strand_profiles = collections.defaultdict(init)

    for match in matches:
        position_info = match.position_info

        relative_offset = int(match.relative_offset)
        strand = position_info[5]
        length = int(position_info[6])

        profile = length_strand_profiles[(length, strand)]

        profile[relative_offset] += 1

    return length_strand_profiles


def get_metagene_profile_df(
    length, type_label, length_strand_profiles, upstream, downstream
):

    reverse_metagene_profile = length_strand_profiles[(length, "-")]
    forward_metagene_profile = length_strand_profiles[(length, "+")]

    metagene_profile = forward_metagene_profile + reverse_metagene_profile[::-1]

    if sum(metagene_profile) == 0:
        return None

    offset = range(-1 * upstream, downstream + 1)

    metagene_profile_df = pd.DataFrame()
    metagene_profile_df["position"] = offset
    metagene_profile_df["count"] = metagene_profile
    metagene_profile_df["type"] = type_label
    metagene_profile_df["length"] = length

    return metagene_profile_df


def main():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""This script extracts the metagene profile
        from reads in a BAM file, possibly filtering by length.""",
    )

    parser.add_argument("bam", help="The bam file")

    parser.add_argument("orfs", help="The annotated transcripts (bed) file")

    parser.add_argument("out", help="The (output) csv.gz counts file")

    parser.add_argument(
        "-p", "--num-cpus", help="The number of processors to use", type=int, default=1
    )

    parser.add_argument(
        "--is-sam",
        help="""If this flag is present, the alignment file will
        be parsed as SAM rather than BAM""",
        action="store_true",
    )

    parser.add_argument(
        "--lengths",
        help="""If specified, then metagene profiles will be
        created for reads of each length. Otherwise, profiles will be created for each
        read length present in the bam file.""",
        type=int,
        nargs="*",
        default=[],
    )

    parser.add_argument(
        "--start-upstream",
        type=int,
        default=default_start_upstream,
        help="""The number of bases upstream of the translation initiation
        site to begin constructing the metagene profile.""",
    )

    parser.add_argument(
        "--start-downstream",
        type=int,
        default=default_start_downstream,
        help="""The number of bases downstream of the translation initiation
        site to end the metagene profile.""",
    )

    parser.add_argument(
        "--end-upstream",
        type=int,
        default=default_end_upstream,
        help="""The number of bases upstream of the translation termination
        site to begin constructing the metagene profile.""",
    )

    parser.add_argument(
        "--end-downstream",
        type=int,
        default=default_end_downstream,
        help="""The number of bases downstream of the translation termination
        site to end the metagene profile.""",
    )

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "[extract-metagene-profiles]: {}".format(" ".join(sys.argv))
    logger.info(msg)

    # first, get the 5' ends of the reads
    alignment_df = bam_utils.get_five_prime_ends(
        args.bam, progress_bar=True, count=True, logger=logger
    )

    msg = "Reading annotations"
    logger.info(msg)
    annotations_df = bed_utils.read_bed(args.orfs, low_memory=False)

    msg = "Constructing canonical translation initiation ORF data frames"
    logger.info(msg)

    m_has_canonical = annotations_df["thick_start"] > -1
    m_forward = annotations_df["strand"] == "+"

    m_canonical_forward = m_has_canonical & m_forward
    m_canonical_reverse = m_has_canonical & ~m_forward

    # forward translation initiation
    start = annotations_df.loc[m_canonical_forward, "thick_start"] - args.start_upstream
    end = annotations_df.loc[m_canonical_forward, "thick_start"] + args.start_downstream
    seqname = annotations_df.loc[m_canonical_forward, "seqname"]
    strand = "+"

    intervals_forward_initiation_bed = get_interval_df(start, end, seqname, strand)

    # reverse translation initation
    start = annotations_df.loc[m_canonical_reverse, "thick_end"] - args.start_downstream
    end = annotations_df.loc[m_canonical_reverse, "thick_end"] + args.start_upstream
    seqname = annotations_df.loc[m_canonical_reverse, "seqname"]
    strand = "-"

    intervals_reverse_initiation_bed = get_interval_df(start, end, seqname, strand)

    # all translation initiation regions
    intervals_initiation_bed = pd.concat(
        [intervals_forward_initiation_bed, intervals_reverse_initiation_bed]
    )

    # make sure we do not double count isoforms with the same starts
    intervals_initiation_bed = intervals_initiation_bed.drop_duplicates()

    msg = "Constructing canonical translation termination ORF data frames"
    logger.info(msg)

    # forward translation termination
    start = annotations_df.loc[m_canonical_forward, "thick_end"] - args.end_upstream
    end = annotations_df.loc[m_canonical_forward, "thick_end"] + args.end_downstream
    seqname = annotations_df.loc[m_canonical_forward, "seqname"]
    strand = "+"

    intervals_forward_termination_bed = get_interval_df(start, end, seqname, strand)

    # reverse translation termination
    start = annotations_df.loc[m_canonical_reverse, "thick_start"] - args.end_downstream
    end = annotations_df.loc[m_canonical_reverse, "thick_start"] + args.end_upstream
    seqname = annotations_df.loc[m_canonical_reverse, "seqname"]
    strand = "-"
    intervals_reverse_termination_bed = get_interval_df(start, end, seqname, strand)

    # all translation termination regions
    intervals_termination_bed = pd.concat(
        [intervals_forward_termination_bed, intervals_reverse_termination_bed]
    )

    # make sure we do not double count isoforms with the same starts
    intervals_termination_bed = intervals_termination_bed.drop_duplicates()

    msg = "Finding translation initiation site matches"
    logger.info(msg)

    initiation_matches = bed_utils.get_all_position_intersections(
        alignment_df, intervals_initiation_bed
    )
    profile_length = args.start_upstream + args.start_downstream + 1
    initiation_length_strand_profiles = get_length_strand_profiles(
        initiation_matches, profile_length
    )

    initiation_keys_str = ",".join(
        str(k) for k in initiation_length_strand_profiles.keys()
    )
    msg = "Initiation keys: {}".format(initiation_keys_str)
    logger.debug(msg)

    msg = "Finding translation termination site matches"
    logger.info(msg)

    termination_matches = bed_utils.get_all_position_intersections(
        alignment_df, intervals_termination_bed
    )
    profile_length = args.end_upstream + args.end_downstream + 1
    termination_length_strand_profiles = get_length_strand_profiles(
        termination_matches, profile_length
    )

    termination_keys_str = ",".join(
        str(k) for k in termination_length_strand_profiles.keys()
    )
    msg = "Termination keys: {}".format(termination_keys_str)
    logger.debug(msg)

    msg = "Extracting metagene profiles"
    logger.info(msg)

    if len(args.lengths) == 0:
        args.lengths = list(alignment_df["length"].unique())

    args.lengths = np.sort(args.lengths)
    args.lengths = [int(l) for l in args.lengths]
    length_str = ",".join(str(l) for l in args.lengths)
    msg = "Profiles will be created for lengths: {}".format(length_str)
    logger.info(msg)

    all_metagene_profile_dfs = []

    for length in args.lengths:
        # first, the profile for this length around initiation sites
        initiation_profile_df = get_metagene_profile_df(
            length,
            "start",
            initiation_length_strand_profiles,
            args.start_upstream,
            args.start_downstream,
        )

        all_metagene_profile_dfs.append(initiation_profile_df)

        # and around termination sites
        termination_profile_df = get_metagene_profile_df(
            length,
            "end",
            termination_length_strand_profiles,
            args.end_upstream,
            args.end_downstream,
        )

        all_metagene_profile_dfs.append(termination_profile_df)

    # filter out all of the profiles which did not have any reads
    all_metagene_profile_dfs = [df for df in all_metagene_profile_dfs if df is not None]

    # join them together in one large data frame
    all_metagene_profile_dfs = pd.concat(all_metagene_profile_dfs)

    msg = "Writing metagene profiles to disk"
    logger.info(msg)
    pandas_utils.write_df(all_metagene_profile_dfs, args.out, index=False)


if __name__ == "__main__":
    main()
