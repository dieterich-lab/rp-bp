#! /usr/bin/env python3

import argparse

import os
import sys
import yaml
import numpy as np
import pandas as pd
import scipy.stats

import pbio.utils.bio as bio
import pbio.utils.bed_utils as bed_utils
import pbio.utils.mygene_utils as mygene_utils
import pbio.misc.math_utils as math_utils
import pbio.misc.parallel as parallel
import pbio.misc.utils as utils
import pbio.misc.pandas_utils as pandas_utils

import pbio.ribo.ribo_utils as ribo_utils
import pbio.ribo.ribo_filenames as filenames

import pyensembl


# TODO: this script causes several SettingWithCopyWarnings. This is expected, so
# we will ignore it for now.
pd.options.mode.chained_assignment = None

import logging
import pbio.misc.logging_utils as logging_utils

logger = logging.getLogger(__name__)

default_fields_to_keep = ["id", "orf_type", "x_1_sum", "x_2_sum", "x_3_sum"]
default_max_micropeptide_len = 300

default_ensembl_release = 79
default_ensembl_species = "mouse"

default_read_filter_percent = 0.25
default_kl_filter_percent = 0.25

default_id_matches = []
default_id_match_names = []
default_overlaps = []
default_overlap_names = []


def add_id_matches(diff_micropeptides, id_file, name):
    msg = "Reading matches id file: {}".format(id_file)
    logger.info(msg)

    matches = pd.read_csv(id_file, header=None, names=["orf_id"])
    matches["orf_id"] = matches["orf_id"].replace("TCONS_", "TCONS")
    matches = set(matches["orf_id"])

    msg = "Finding matches"
    logger.info(msg)
    m_match_a = diff_micropeptides["A"].isin(matches)
    m_match_b = diff_micropeptides["B"].isin(matches)

    match_name_a = "{}_A".format(name)
    match_name_b = "{}_B".format(name)

    diff_micropeptides[match_name_a] = "No"
    diff_micropeptides[match_name_b] = "No"

    diff_micropeptides.loc[m_match_a, match_name_a] = "Yes"
    diff_micropeptides.loc[m_match_b, match_name_b] = "Yes"

    return diff_micropeptides


def add_overlaps(diff_micropeptides, overlap_file, name, bed_df_a, bed_df_b, exons):
    msg = "Reading overlaps file: {}".format(overlap_file)
    logger.info(msg)
    overlap_bed = bed_utils.read_bed(overlap_file)

    msg = "Finding overlaps"
    a_overlaps = bed_utils.get_bed_overlaps(bed_df_a, overlap_bed, exons=exons)
    a_overlaps_ids = {to.a_info for to in a_overlaps}

    b_overlaps = bed_utils.get_bed_overlaps(bed_df_b, overlap_bed, exons=exons)
    b_overlaps_ids = {to.a_info for to in b_overlaps}

    m_match_a = diff_micropeptides["A"].isin(a_overlaps_ids)
    m_match_b = diff_micropeptides["B"].isin(b_overlaps_ids)

    match_name_a = "{}_A".format(name)
    match_name_b = "{}_B".format(name)

    diff_micropeptides[match_name_a] = "No"
    diff_micropeptides[match_name_b] = "No"

    diff_micropeptides.loc[m_match_a, match_name_a] = "Yes"
    diff_micropeptides.loc[m_match_b, match_name_b] = "Yes"

    return diff_micropeptides


X_I_FIELDS = ["x_1_sum", "x_2_sum", "x_3_sum"]


def get_orf_kl_entry(overlap, overlap_type, bf_df_a, bf_df_b):

    # grab the right ORFs
    m_a_id = bf_df_a["id"] == overlap.a_info
    m_b_id = bf_df_b["id"] == overlap.b_info

    # pull out the frame counts
    a_counts = bf_df_a.loc[m_a_id, X_I_FIELDS].iloc[0]
    b_counts = bf_df_b.loc[m_b_id, X_I_FIELDS].iloc[0]

    a_counts = np.array(a_counts)
    b_counts = np.array(b_counts)

    # laplacian smoothing
    a_counts += 1
    b_counts += 1

    kl = math_utils.calculate_symmetric_kl_divergence(
        a_counts, b_counts, scipy.stats.entropy
    )

    ret = {
        "A": overlap.a_info,
        "B": overlap.b_info,
        "kl": kl,
        "overlap_type": overlap_type,
    }

    return ret


def get_overlap_df(a, b, overlap_type, bf_df_a, bf_df_b):
    # find overlapping micropeptides
    overlap = bed_utils.get_bed_overlaps(a, b)

    overlap_df = parallel.apply_iter_simple(
        overlap, get_orf_kl_entry, overlap_type, bf_df_a, bf_df_b, progress_bar=True
    )

    overlap_df = pd.DataFrame(overlap_df)

    return overlap_df


def get_transcript_and_biotype(row, condition, ensembl, ensembl_transcript_ids):
    mackowiak_id = row[condition]
    if mackowiak_id is None:
        return None

    info = bio.parse_mackowiak_id(mackowiak_id)
    transcript_id, seqname, start, end, strand = info

    transcript_id = transcript_id
    biotype = None
    gene_id = None

    if transcript_id in ensembl.transcript_ids():
        t = ensembl.transcript_by_id(transcript_id)
        biotype = t.biotype
        gene_id = t.gene_id
    else:
        # search based on the location
        t = ensembl.transcript_ids_at_locus(seqname, start, end, strand)

        if len(t) > 0:
            transcript_id = t[0]
            t = ensembl.transcript_by_id(t[0])
            biotype = t.biotype
            gene_id = t.gene_id

    bt = "biotype_{}".format(condition)
    ti = "transcript_id_{}".format(condition)
    gi = "gene_id_{}".format(condition)

    ret = {condition: mackowiak_id, bt: biotype, ti: transcript_id, gi: gene_id}

    return ret


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script extracts the differential micropeptides from two "
        "conditions. Please see the documentation in redmine for more details.\n\n"
        "Please see the pyensembl (https://github.com/hammerlab/pyensembl) "
        "documentation for more information about the ensembl release and species.",
    )

    parser.add_argument("config", help="The (yaml) config file")
    parser.add_argument("name_a", help="The name of the first condition")
    parser.add_argument("name_b", help="The name of the second condition")
    parser.add_argument("out", help="The output (.csv.gz or .xlsx) file")

    parser.add_argument(
        "-a",
        "--append-sheet",
        help="If this flag is given, "
        "then a worksheet with the name '<name_a>,<name_b>' will be appended "
        "to the .xlsx file given by out (if it exists)",
        action="store_true",
    )

    parser.add_argument(
        "-f",
        "--filter",
        help="If this flag is present, then "
        "the output will be filtered to include only the differential "
        "micropeptides with the highest KL-divergence and read coverage",
        action="store_true",
    )

    parser.add_argument(
        "--read-filter-percent",
        help="If the the --filter flag "
        "is given, then only the top --read-filter-percent micropeptides will "
        "be considered for the final output. They still must meet the KL-"
        "divergence filtering criteria.",
        type=float,
        default=default_read_filter_percent,
    )

    parser.add_argument(
        "--kl-filter-percent",
        help="If the the --filter flag "
        "is given, then only the top --read-kl-percent micropeptides will "
        "be considered for the final output. They still must meet the read "
        "coverage filtering criteria.",
        type=float,
        default=default_kl_filter_percent,
    )

    parser.add_argument(
        "--id-matches",
        help="This is a list of files which "
        "contain ORF identifiers to compare to the differential micropeptides. "
        "For each of the files given, two columns will be added to the output "
        "which indicate if either A or B appear in the respective file. Each "
        "file should have a single ORF identifier on each line and contain "
        "nothing else.",
        nargs="*",
        default=default_id_matches,
    )

    parser.add_argument(
        "--id-match-names",
        help="A name to include in the "
        "output file for each --id-matches file. The number of names must "
        "match the number of files.",
        nargs="*",
        default=default_id_match_names,
    )

    parser.add_argument(
        "--overlaps",
        help="This is a list of bed12+ files "
        "which will be compared to the differential micropeptides. Two columns "
        "(one for A, one for B) will be added to the output which indicate if "
        "the respective micropeptides overlap a feature in each file by at "
        "least 1 bp.",
        nargs="*",
        default=default_overlaps,
    )

    parser.add_argument(
        "--overlap-names",
        help="A name to include in the "
        "output file for each --overlaps file. The number of names must match "
        "the number of files.",
        nargs="*",
        default=default_overlap_names,
    )

    parser.add_argument(
        "-r",
        "--ensembl-release",
        help="The version of Ensembl "
        "to use when mapping transcript identifiers to gene identifiers",
        type=int,
        default=default_ensembl_release,
    )

    parser.add_argument(
        "-s",
        "--ensembl-species",
        help="The Ensembl species "
        "to use when mapping transcript identifiers to gene identifiers",
        default=default_ensembl_species,
    )

    parser.add_argument(
        "--a-is-single-sample",
        help="By default, this script "
        "assumes the predictions come from merged replicates. If name_a is from "
        "a single sample, this flag should be given. It is necessary to find "
        "the correct filenames.",
        action="store_true",
    )

    parser.add_argument(
        "--b-is-single-sample",
        help="By default, this script "
        "assumes the predictions come from merged replicates. If name_b is from "
        "a single sample, this flag should be given. It is necessary to find "
        "the correct filenames.",
        action="store_true",
    )

    parser.add_argument(
        "--fields-to-keep",
        help="The fields to keep from the " "Bayes factor file for each condition",
        nargs="*",
        default=default_fields_to_keep,
    )

    parser.add_argument(
        "--max-micropeptide-len",
        help="The maximum (inclusive) " "length of ORFs considered as micropeptides",
        type=int,
        default=default_max_micropeptide_len,
    )

    parser.add_argument(
        "--do-not-fix-tcons",
        help='By default, the "TCONS_" '
        "identifiers from StringTie, etc., do not parse correctly; this script "
        "update the identifiers so that will parse correctly unless instructed not "
        "to. The script is likely to crash if the identifiers are not fixed.",
        action="store_true",
    )

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "Loading ensembl database"
    logger.info(msg)

    ensembl = pyensembl.EnsemblRelease(
        release=args.ensembl_release, species=args.ensembl_species
    )
    ensembl.db

    msg = "Checking the id-match and overlaps files"
    logger.info(msg)

    if len(args.id_matches) != len(args.id_match_names):
        msg = (
            "The number of --id-matches files and --id-match-names do not "
            "match. {} files and {} names".format(
                len(args.id_matches), len(args.id_match_names)
            )
        )

        raise ValueError(msg)

    if len(args.overlaps) != len(args.overlap_names):
        msg = (
            "The number of --overlaps files and --overlaps-names do not "
            "match. {} files and {} names".format(
                len(args.overlaps), len(args.overlap_names)
            )
        )

        raise ValueError(msg)

    utils.check_files_exist(args.id_matches)
    utils.check_files_exist(args.overlaps)

    if args.filter:
        msg = "Validating filter percentages"
        logger.info(msg)

        math_utils.check_range(
            args.read_filter_percent, 0, 1, variable_name="--read-filter-percent"
        )

        math_utils.check_range(
            args.kl_filter_percent, 0, 1, variable_name="--kl-filter-percent"
        )

    msg = "Extracting file names"
    logger.info(msg)

    config = yaml.load(open(args.config), Loader=yaml.FullLoader)

    note_str = config.get("note", None)

    # keep multimappers?
    is_unique = not ("keep_riboseq_multimappers" in config)

    # and the smoothing parameters
    fraction = config.get("smoothing_fraction", None)
    reweighting_iterations = config.get("smoothing_reweighting_iterations", None)

    lengths_a = None
    offsets_a = None

    if args.a_is_single_sample:
        lengths_a, offsets_a = ribo_utils.get_periodic_lengths_and_offsets(
            config, args.name_a, is_unique=is_unique
        )

    bayes_factors_a = filenames.get_riboseq_bayes_factors(
        config["riboseq_data"],
        args.name_a,
        length=lengths_a,
        offset=offsets_a,
        is_unique=is_unique,
        note=note_str,
        fraction=fraction,
        reweighting_iterations=reweighting_iterations,
    )

    if not os.path.exists(bayes_factors_a):
        msg = "Could not find the Bayes factor file for {}. ({}). Quitting.".format(
            args.name_a, bayes_factors_a
        )
        raise FileNotFoundError(msg)

    predicted_orfs_a = filenames.get_riboseq_predicted_orfs(
        config["riboseq_data"],
        args.name_a,
        length=lengths_a,
        offset=offsets_a,
        is_unique=is_unique,
        note=note_str,
        fraction=fraction,
        reweighting_iterations=reweighting_iterations,
        is_filtered=True,
        is_chisq=False,
    )

    if not os.path.exists(predicted_orfs_a):
        msg = "Could not find the predictions bed file for {}. ({}). Quitting.".format(
            args.name_a, predicted_orfs_a
        )
        raise FileNotFoundError(msg)

    lengths_b = None
    offsets_b = None
    if args.b_is_single_sample:
        lengths_b, offsets_b = ribo_utils.get_periodic_lengths_and_offsets(
            config, args.name_b, is_unique=is_unique
        )

    bayes_factors_b = filenames.get_riboseq_bayes_factors(
        config["riboseq_data"],
        args.name_b,
        length=lengths_b,
        offset=offsets_b,
        is_unique=is_unique,
        note=note_str,
        fraction=fraction,
        reweighting_iterations=reweighting_iterations,
    )

    if not os.path.exists(bayes_factors_b):
        msg = "Could not find the Bayes factor file for {}. ({}). Quitting.".format(
            args.name_b, bayes_factors_b
        )
        raise FileNotFoundError(msg)

    predicted_orfs_b = filenames.get_riboseq_predicted_orfs(
        config["riboseq_data"],
        args.name_b,
        length=lengths_b,
        offset=offsets_b,
        is_unique=is_unique,
        note=note_str,
        fraction=fraction,
        reweighting_iterations=reweighting_iterations,
        is_filtered=True,
        is_chisq=False,
    )

    if not os.path.exists(predicted_orfs_b):
        msg = "Could not find the predictions bed file for {}. ({}). Quitting.".format(
            args.name_b, predicted_orfs_b
        )
        raise FileNotFoundError(msg)

    exons_file = filenames.get_exons(
        config["genome_base_path"], config["genome_name"], note=config.get("orf_note")
    )

    if not os.path.exists(exons_file):
        msg = "Could not find the exons file ({}). Quitting.".format(exons_file)
        raise FileNotFoundError(msg)

    msg = "Reading the exons"
    logger.info(msg)

    exons = bed_utils.read_bed(exons_file)

    msg = "Reading the BF files"
    logger.info(msg)

    bf_df_a = bed_utils.read_bed(bayes_factors_a)
    bf_df_b = bed_utils.read_bed(bayes_factors_b)

    msg = "Reading the predictions files"
    logger.info(msg)

    bed_df_a = bed_utils.read_bed(predicted_orfs_a)
    bed_df_b = bed_utils.read_bed(predicted_orfs_b)

    differential_micropeptide_dfs = []

    # extract micropeptides
    msg = "Extracting micropeptides"
    logger.info(msg)

    m_micropeptides_a = bed_df_a["orf_len"] <= args.max_micropeptide_len
    m_micropeptides_b = bed_df_b["orf_len"] <= args.max_micropeptide_len

    micropeptides_a = bed_df_a[m_micropeptides_a]
    micropeptides_b = bed_df_b[m_micropeptides_b]

    long_orfs_a = bed_df_a[~m_micropeptides_a]
    long_orfs_b = bed_df_b[~m_micropeptides_b]

    msg = "Finding micropeptides in A with no overlap in B"
    logger.info(msg)

    micropeptides_a_no_match_b = bed_utils.subtract_bed(
        micropeptides_a, bed_df_b, exons=exons
    )

    micropeptides_a_no_match_b_df = pd.DataFrame()
    micropeptides_a_no_match_b_df["A"] = list(micropeptides_a_no_match_b)
    micropeptides_a_no_match_b_df["B"] = None
    micropeptides_a_no_match_b_df["kl"] = np.inf
    micropeptides_a_no_match_b_df["overlap_type"] = "micro_a_only"

    differential_micropeptide_dfs.append(micropeptides_a_no_match_b_df)

    msg = "Finding micropeptides in B with no overlap in A"
    logger.info(msg)

    micropeptides_b_no_match_a = bed_utils.subtract_bed(
        micropeptides_b, bed_df_a, exons=exons
    )

    micropeptides_b_no_match_a_df = pd.DataFrame()
    micropeptides_b_no_match_a_df["B"] = list(micropeptides_b_no_match_a)
    micropeptides_b_no_match_a_df["A"] = None
    micropeptides_b_no_match_a_df["kl"] = np.inf
    micropeptides_b_no_match_a_df["overlap_type"] = "micro_b_only"

    differential_micropeptide_dfs.append(micropeptides_b_no_match_a_df)

    msg = "Finding overlapping micropeptides"
    logger.info(msg)

    micropeptides_a_micropeptides_b_df = get_overlap_df(
        micropeptides_a, micropeptides_b, "micro_a_micro_b", bf_df_a, bf_df_b
    )
    differential_micropeptide_dfs.append(micropeptides_a_micropeptides_b_df)

    micropeptides_a_long_b_df = get_overlap_df(
        micropeptides_a, long_orfs_b, "micro_a_long_b", bf_df_a, bf_df_b
    )
    differential_micropeptide_dfs.append(micropeptides_a_long_b_df)

    micropeptides_b_long_a_df = get_overlap_df(
        long_orfs_a, micropeptides_b, "long_a_micro_b", bf_df_a, bf_df_b
    )
    differential_micropeptide_dfs.append(micropeptides_b_long_a_df)

    differential_micropeptides_df = pd.concat(differential_micropeptide_dfs)

    msg = "Adding read count information"
    logger.info(msg)

    res = differential_micropeptides_df.merge(
        bf_df_a[args.fields_to_keep], left_on="A", right_on="id", how="left"
    )
    to_rename = {f: "{}_A".format(f) for f in args.fields_to_keep}
    res = res.rename(columns=to_rename)
    res = res.drop("id_A", axis=1)

    res = res.merge(
        bf_df_b[args.fields_to_keep], left_on="B", right_on="id", how="left"
    )
    to_rename = {f: "{}_B".format(f) for f in args.fields_to_keep}
    res = res.rename(columns=to_rename)
    res = res.drop("id_B", axis=1)

    id_columns = ["A", "B"]
    res = res.drop_duplicates(subset=id_columns)

    if not args.do_not_fix_tcons:
        # replace TCONS_ with TCONS
        res["A"] = res["A"].str.replace("TCONS_", "TCONS")
        res["B"] = res["B"].str.replace("TCONS_", "TCONS")

    msg = "Extracting the genes and their biotypes using pyensembl"
    logger.info(msg)

    ensembl = pyensembl.EnsemblRelease(
        release=args.ensembl_release, species=args.ensembl_species
    )
    ensembl_transcript_ids = set(ensembl.transcript_ids())

    biotypes_a = parallel.apply_df_simple(
        res, get_transcript_and_biotype, "A", ensembl, ensembl_transcript_ids
    )
    biotypes_b = parallel.apply_df_simple(
        res, get_transcript_and_biotype, "B", ensembl, ensembl_transcript_ids
    )

    biotypes_a = utils.remove_nones(biotypes_a)
    biotypes_b = utils.remove_nones(biotypes_b)

    biotypes_a = pd.DataFrame(biotypes_a)
    biotypes_b = pd.DataFrame(biotypes_b)

    res = res.merge(biotypes_a, on="A", how="left")
    res = res.merge(biotypes_b, on="B", how="left")

    msg = "Pulling annotations from mygene.info"
    logger.info(msg)

    # pull annotations from mygene
    gene_info_a = mygene_utils.query_mygene(res["gene_id_A"])
    gene_info_b = mygene_utils.query_mygene(res["gene_id_B"])

    # and add the mygene info
    res = res.merge(gene_info_a, left_on="gene_id_A", right_on="gene_id", how="left")

    to_rename = {f: "{}_A".format(f) for f in gene_info_a.columns}
    to_rename.pop("gene_id")
    res = res.rename(columns=to_rename)
    res = res.drop("gene_id", axis=1)

    res = res.merge(gene_info_b, left_on="gene_id_B", right_on="gene_id", how="left")

    to_rename = {f: "{}_B".format(f) for f in gene_info_a.columns}
    to_rename.pop("gene_id")
    res = res.rename(columns=to_rename)
    res = res.drop("gene_id", axis=1)

    msg = "Removing duplicates"
    logger.info(msg)
    id_columns = ["A", "B"]
    res = res.drop_duplicates(subset=id_columns)

    msg = "Adding --id-matches columns"
    logger.info(msg)

    for (id_match_file, name) in zip(args.id_matches, args.id_match_names):
        res = add_id_matches(res, id_match_file, name)

    msg = "Adding --overlaps columns"
    logger.info(msg)

    for (overlap_file, name) in zip(args.overlaps, args.overlap_names):
        res = add_overlaps(res, overlap_file, name, bed_df_a, bed_df_b, exons)

    msg = "Sorting by in-frame reads"
    logger.info(msg)

    res["x_1_sum_A"] = res["x_1_sum_A"].fillna(0)
    res["x_1_sum_B"] = res["x_1_sum_B"].fillna(0)
    res["x_1_sum"] = res["x_1_sum_A"] + res["x_1_sum_B"]
    res = res.sort_values("x_1_sum", ascending=False)

    if args.filter:
        msg = "Filtering the micropeptides by read coverage and KL-divergence"
        logger.info(msg)

        x_1_sum_ranks = res["x_1_sum"].rank(
            method="min", na_option="top", ascending=False
        )
        num_x_1_sum_ranks = x_1_sum_ranks.max()
        max_good_x_1_sum_rank = num_x_1_sum_ranks * args.read_filter_percent
        m_good_x_1_sum_rank = x_1_sum_ranks <= max_good_x_1_sum_rank

        msg = "Number of micropeptides passing read filter: {}".format(
            sum(m_good_x_1_sum_rank)
        )
        logger.debug(msg)

        kl_ranks = res["kl"].rank(method="dense", na_option="top", ascending=False)
        num_kl_ranks = kl_ranks.max()
        max_good_kl_rank = num_kl_ranks * args.kl_filter_percent
        m_good_kl_rank = kl_ranks <= max_good_kl_rank

        msg = "Number of micropeptides passing KL filter: {}".format(
            sum(m_good_kl_rank)
        )
        logger.debug(msg)

        m_both_filters = m_good_x_1_sum_rank & m_good_kl_rank

        msg = "Number of micropeptides passing both filters: {}".format(
            sum(m_both_filters)
        )
        logger.debug(msg)

        res = res[m_both_filters]

    msg = "Writing differential micropeptides to disk"
    logger.info(msg)

    if args.append_sheet is None:
        pandas_utils.write_df(res, args.out, index=False)
    else:
        sheet_name = "{},{}".format(args.name_a, args.name_b)
        utils.append_to_xlsx(res, args.out, sheet=sheet_name, index=False)


if __name__ == "__main__":
    main()
