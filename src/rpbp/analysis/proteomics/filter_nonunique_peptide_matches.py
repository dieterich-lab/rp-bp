#! /usr/bin/env python3

import argparse
import pandas as pd

import pbiotools.misc.utils as utils
import pbiotools.misc.pandas_utils as pandas_utils
import pbiotools.misc.parallel as parallel

import logging
import pbiotools.misc.logging_utils as logging_utils

logger = logging.getLogger(__name__)

default_num_cpus = 1


def parse_matches(row):
    """This helper function splits the peptide matches into a list of
    individual matches.
    """
    if row["num_matches"] == 0:
        return None

    peptide_matches = row["peptide_matches"].split(";")

    ret = [
        {
            "orf_id": row["orf_id"],
            "orf_sequence": row["orf_sequence"],
            "peptide": peptide,
        }
        for peptide in peptide_matches
    ]

    return ret


def merge_group(g):
    """This helper function merges the individual (unique) matches back
    into a single record.
    """
    orf_id = g.iloc[0]["orf_id"]
    orf_sequence = g.iloc[0]["orf_sequence"]
    num_matches = len(g)
    peptide_matches = ";".join(g["peptide"])

    ret = {
        "orf_id": orf_id,
        "orf_sequence": orf_sequence,
        "num_matches": num_matches,
        "peptide_matches": peptide_matches,
    }

    return ret


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script removes all of the peptides which match to multiple "
        "ORFs from the results found with get-all-orf-peptide-matches.",
    )

    parser.add_argument(
        "peptide_matches",
        help="The peptide matches file produced " "by get-all-orf-peptide-matches",
    )
    parser.add_argument(
        "out",
        help="A similar peptide matches file which "
        "contains only peptides which match to a unique ORF",
    )

    parser.add_argument(
        "-p",
        "--num-cpus",
        help="The number of CPUs to use",
        type=int,
        default=default_num_cpus,
    )

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "Reading peptide matches"
    logger.info(msg)

    peptide_matches = pd.read_csv(args.peptide_matches)

    msg = "Splitting the grouped matches into individual peptide matches"
    logger.info(msg)

    matches = parallel.apply_parallel(
        peptide_matches, args.num_cpus, parse_matches, progress_bar=True
    )

    msg = "Removing peptides which match to multiple ORFs"
    logger.info(msg)

    matches = utils.remove_nones(matches)
    matches = utils.flatten_lists(matches)
    matches_df = pd.DataFrame(matches)
    unique_matches_df = matches_df.drop_duplicates(subset="peptide", keep=False)

    msg = "Merging the ORF-peptide matches back to single records"
    logger.info(msg)

    unique_groups = unique_matches_df.groupby("orf_id")
    merged_unique_groups = parallel.apply_parallel_groups(
        unique_groups, args.num_cpus, merge_group, progress_bar=True
    )

    merged_unique_df = pd.DataFrame(merged_unique_groups)

    msg = "Re-adding the ORFs which no longer have peptide matches"
    logger.info(msg)

    m_still_has_match = peptide_matches["orf_id"].isin(merged_unique_df["orf_id"])
    peptide_matches.loc[~m_still_has_match, "num_matches"] = 0
    peptide_matches.loc[~m_still_has_match, "peptide_matches"] = 0

    peps = [merged_unique_df, peptide_matches[~m_still_has_match]]
    merged_unique_df = pd.concat(peps)

    msg = "Writing the ORFs with unique matches to disk"
    logger.info(msg)

    pandas_utils.write_df(merged_unique_df, args.out, index=False)


if __name__ == "__main__":
    main()
