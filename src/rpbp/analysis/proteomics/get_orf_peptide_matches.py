#! /usr/bin/env python3

import argparse
import pandas as pd
import logging

import sys

import pbiotools.utils.fastx_utils as fastx_utils
import pbiotools.misc.logging_utils as logging_utils
import pbiotools.misc.parallel as parallel
import pbiotools.misc.pandas_utils as pandas_utils

logger = logging.getLogger(__name__)

default_num_cpus = 1

default_num_groups = 100

default_num_peptides = 0

default_peptide_separator = "\t"
default_peptide_filter_field = "PEP"
default_peptide_filter_value = 0.1


def get_match_series(o, peptide):
    """This function extracts relevant information from the orf and peptide
    string. It returns the results as a pd.Series.
    """
    ret = {"peptide": peptide, "orf_id": o["orf_id"]}
    return pd.Series(ret)


def find_matching_orfs(peptide, orfs):
    """This function finds all of the ORFs which include the peptide as an
    EXACT SUBSTRING.

    Args:
        peptide (pd.Series): the peptide to search for

        orfs (pd.DataFrame): All of the predicted ORFs in which to search

    Returns:
        pd.DataFrame: containing the peptide, orf_id and orf_sequence
            of all matches.
    """
    peptide_seq = peptide["Sequence"]
    mask_matching = orfs["orf_sequence"].str.contains(peptide_seq)

    # short-circuit, when possible
    if sum(mask_matching) == 0:
        return None

    matching_orfs = orfs[mask_matching]
    ret = matching_orfs.apply(get_match_series, args=(peptide_seq,), axis=1)
    return ret


def find_matching_orfs_group(peptides, orfs):
    """A helper function to call find_matching_orfs on a pd.GroupBy of peptides."""
    ret = parallel.apply_df_simple(peptides, find_matching_orfs, orfs)
    # progress_bar=True)
    ret = [r for r in ret if r is not None]
    if len(ret) == 0:
        return None
    return pd.concat(ret)


def count_matches(peptide_matches):
    """This function counts the number of matches in the given group. It returns
    a series containing the orf_id and number of matches.
    """
    num_matches = len(peptide_matches)
    peptide_matches_str = ";".join(peptide_matches["peptide"])
    orf_id = peptide_matches.iloc[0]["orf_id"]
    ret = {
        "num_matches": num_matches,
        "peptide_matches": peptide_matches_str,
        "orf_id": orf_id,
    }
    return pd.Series(ret)


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script uses the peptides.txt file from MaxQuant to determine "
        "which predicted ORFs have some proteomics evidence.\n\nIt contains "
        "some hard-coded field names.",
    )
    parser.add_argument(
        "predicted_proteins", help="The (fasta, protein) file of " "predicted ORFs"
    )
    parser.add_argument("peptides", help="The peptides.txt file produced by MaxQuant")
    parser.add_argument(
        "out",
        help="The output (csv.gz) file containing the predicted "
        "ORFs and their coverage",
    )

    parser.add_argument(
        "--num-cpus",
        help="The number of CPUs to use for searching",
        type=int,
        default=default_num_cpus,
    )

    parser.add_argument(
        "--peptide-filter-field",
        help="The field to use for filtering " "the peptides from MaxQuant",
        default=default_peptide_filter_field,
    )
    parser.add_argument(
        "--peptide-filter-value",
        help="All peptides with a value greater "
        "than the filter value will be removed",
        type=float,
        default=default_peptide_filter_value,
    )

    parser.add_argument(
        "--peptide-separator",
        help="The separator in the --peptide file",
        default=default_peptide_separator,
    )

    parser.add_argument(
        "-g",
        "--num-groups",
        help="The number of groups into which to split "
        "the ORFs. More groups means the progress bar is updated more frequently but incurs "
        "more overhead because of the parallel calls.",
        type=int,
        default=default_num_groups,
    )

    parser.add_argument(
        "--num-peptides",
        help="If n>0, then only the first n peptide "
        "sequences will be used to calculate coverage. This is for testing.",
        type=int,
        default=default_num_peptides,
    )

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "[get-orf-peptide-matches]: {}".format(" ".join(sys.argv))
    logger.info(msg)

    msg = "Reading and filtering peptides"
    logger.info(msg)

    peptides = pd.read_csv(args.peptides, sep=args.peptide_separator)
    mask_filter = peptides[args.peptide_filter_field] < args.peptide_filter_value
    peptides = peptides[mask_filter]
    peptide_sequences = pd.DataFrame(peptides["Sequence"])

    if args.num_peptides > 0:
        peptide_sequences = peptide_sequences.head(args.num_peptides)

    msg = "Number of filtered peptides: {}".format(len(peptide_sequences))
    logger.info(msg)

    msg = "Reading predicted ORFs into a data frame"
    logger.info(msg)

    # TODO: use read iterator
    predicted_orfs = fastx_utils.get_read_iterator(args.predicted_proteins)
    orf_ids = []
    orf_sequences = []

    for orf_id, seq in predicted_orfs:
        orf_ids.append(orf_id)
        orf_sequences.append(seq)

    predicted_orfs_df = pd.DataFrame()
    predicted_orfs_df["orf_id"] = orf_ids
    predicted_orfs_df["orf_sequence"] = orf_sequences

    msg = "Searching for matching peptides"
    logger.info(msg)

    peptide_matches = parallel.apply_parallel_split(
        peptide_sequences,
        args.num_cpus,
        find_matching_orfs_group,
        predicted_orfs_df,
        progress_bar=True,
        num_groups=args.num_groups,
    )

    # filter out the Nones to avoid DataFrame conversion problems
    msg = "Joining results back into large data frame"
    logger.info(msg)

    peptide_matches = [pm for pm in peptide_matches if pm is not None]
    peptide_matches = pd.concat(peptide_matches)

    # now, we have a data frame of matches (fields: peptide, orf_id)
    msg = "Getting peptide coverage of ORFs"
    logger.info(msg)

    # first, count the matches for each ORF
    peptide_matches_groups = peptide_matches.groupby("orf_id")

    orf_matches = parallel.apply_parallel_groups(
        peptide_matches_groups, args.num_cpus, count_matches, progress_bar=True
    )
    orf_matches = pd.DataFrame(orf_matches)

    # then join back on the original list of ORFs to have entries for ORFs
    # with no peptide matches
    predicted_orf_coverage = pd.merge(
        predicted_orfs_df, orf_matches, on="orf_id", how="left"
    )

    # and patch the holes in the data frame
    predicted_orf_coverage = predicted_orf_coverage.fillna(0)

    msg = "Writing coverage information to disk"
    pandas_utils.write_df(predicted_orf_coverage, args.out, index=False)


if __name__ == "__main__":
    main()
