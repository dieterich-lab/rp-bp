#! /usr/bin/env python3

import argparse
import functools
import gc
import logging
import sys
import numpy as np
import scipy.sparse
import tqdm

import pbiotools.utils.bed_utils as bed_utils
import pbiotools.misc.logging_utils as logging_utils
import pbiotools.misc.pandas_utils as pandas_utils

import pbiotools.misc.math_utils as math_utils
import pbiotools.misc.utils as utils
import pbiotools.misc.parallel as parallel
import pbiotools.misc.slurm as slurm

import rpbp.ribo_utils.utils as ribo_utils

from rpbp.defaults import default_num_groups

logger = logging.getLogger(__name__)

# --num-exons is not used in the the Rp-Bp pipeline
# --seqname-prefix can be overridden via config arguments, but not
# specified in defaults, as mostly unused


def get_p_site_intersections(seqname, strand, p_sites, exons_df):

    # only the things in the right direction, etc.
    m_exons_seqname = exons_df["seqname"] == seqname
    m_p_sites_seqname = p_sites["seqname"] == seqname

    m_exons_strand = exons_df["strand"] == strand
    m_p_sites_strand = p_sites["strand"] == strand

    p_site_positions = p_sites.loc[m_p_sites_seqname & m_p_sites_strand, "start"]
    exon_starts = exons_df.loc[m_exons_seqname & m_exons_strand, "start"]
    exon_ends = exons_df.loc[m_exons_seqname & m_exons_strand, "end"]

    p_site_positions = np.array(p_site_positions)
    exon_starts = np.array(exon_starts)
    exon_ends = np.array(exon_ends)

    exon_info_fields = ["transcript_start", "orf_num"]
    exon_info = np.array(
        exons_df.loc[m_exons_seqname & m_exons_strand, exon_info_fields]
    )

    intersections = bed_utils.get_position_intersections(
        p_site_positions, exon_starts, exon_ends, exon_info
    )

    return intersections


def get_all_p_site_intersections(exons_psites, num_orfs, max_orf_len):
    """This function finds the intersection of p_sites across all seqnames
    and strands in the given set of exons. It creates a sparse matrix
    representing the corresponding ORF profiles and returns that.
    """

    exons_df, p_sites = exons_psites
    seqnames = exons_df["seqname"].unique()
    strands = ("+", "-")

    msg = "Unique sequences: {}".format(seqnames)
    logger.debug(msg)

    profiles = scipy.sparse.lil_matrix((num_orfs, max_orf_len))

    for seqname in seqnames:
        for strand in strands:
            p_site_intersections = get_p_site_intersections(
                seqname, strand, p_sites, exons_df
            )

            for intersection in p_site_intersections:
                orf_num = intersection[2][1]

                p_site_pos = intersection[1] + intersection[2][0]

                # object array values taken as numpy.float64 ...
                profiles[orf_num, p_site_pos.astype(int)] += 1

    return profiles.tocsr()


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""This script constructs the profile for each ORF. It
        first adjusts the mapped read positions to properly align with the P-sites. Second, it uses
        a custom chrom-sweep algorithm to find the coverage of each position in each exon of each ORF.
        Finally, the ORF exons are glued together to find the profile of the entire ORF.""",
    )

    parser.add_argument(
        "bam", help="The bam file including filtered (unique, etc.) alignments"
    )

    parser.add_argument("orfs", help="The (bed12) file containing the ORFs")

    parser.add_argument("exons", help="The (bed6+2) file containing the exons")

    parser.add_argument(
        "out", help="The (mtx.gz) output file containing the ORF profiles"
    )

    parser.add_argument(
        "-l",
        "--lengths",
        help="""If any values are given, then only reads which have
        those lengths will be included in the signal construction.""",
        type=int,
        default=[],
        nargs="*",
    )

    parser.add_argument(
        "-o",
        "--offsets",
        help="""The 5' end of reads will be shifted by this amount.
        There must be one offset value for each length (given by the --lengths argument.""",
        type=int,
        default=[],
        nargs="*",
    )

    parser.add_argument(
        "-k",
        "--num-exons",
        help="If k>0, then only the first k exons will be processed.",
        type=int,
        default=0,
    )

    parser.add_argument(
        "-g",
        "--num-groups",
        help=""""The number of groups into which to split the exons.
        More groups means the progress bar is updated more frequently but incurs more overhead because of the
        parallel calls.""",
        type=int,
        default=default_num_groups,
    )

    parser.add_argument(
        "--seqname-prefix",
        help="""If present, this string will be prepended to the
        seqname field of the ORFs.""",
        default="",
    )

    slurm.add_sbatch_options(parser)
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "[extract-orf-profiles]: {}".format(" ".join(sys.argv))
    logger.info(msg)

    # make sure the number of lengths and offsets match
    if len(args.lengths) != len(args.offsets):
        msg = "The number of --lengths and --offsets do not match."
        raise ValueError(msg)

    # make sure the necessary files exist
    required_files = [args.bam, args.orfs, args.exons]
    msg = "[extract-orf-profiles]: Some input files were missing: "
    utils.check_files_exist(required_files, msg=msg)

    if args.use_slurm:
        cmd = " ".join(sys.argv)
        slurm.check_sbatch(cmd, args=args)
        return

    msg = "Finding P-sites"
    logger.info(msg)

    p_sites = ribo_utils.get_p_sites(args.bam, args.lengths, args.offsets)

    # we do not need the data frame anymore, so save some memory
    msg = "Reading exons"
    logger.info(msg)
    exons = bed_utils.read_bed(args.exons, low_memory=False)

    msg = "Reading ORFs"
    logger.info(msg)

    orfs = bed_utils.read_bed(args.orfs, low_memory=False)

    if len(args.seqname_prefix) > 0:
        orfs["seqname"] = args.seqname_prefix + orfs["seqname"]
        exons["seqname"] = args.seqname_prefix + exons["seqname"]

    if args.num_exons > 0:
        exons = exons.head(args.num_exons)

    num_orfs = orfs["orf_num"].max() + 1
    max_orf_len = orfs["orf_len"].max()

    msg = "Adding the ORF index to the exons"
    logger.info(msg)

    orf_fields = ["id", "orf_num"]
    exons_orfs = exons.merge(orfs[orf_fields], on="id")

    msg = "Splitting exons and P-sites"
    logger.info(msg)
    exon_groups = pandas_utils.split_df(exons_orfs, args.num_groups)

    exons_dfs = []
    psites_dfs = []

    for group_index, exon_group in exon_groups:
        # pull out only the p-sites that come from these chromosomes
        seqnames = set(exon_group["seqname"].unique())
        m_psites = p_sites["seqname"].isin(seqnames)

        exons_dfs.append(exon_group)
        psites_dfs.append(p_sites[m_psites])

    # we no longer need the full list of psites
    del p_sites
    del exons_orfs
    del exon_groups
    del exons
    gc.collect()
    exons_psites = zip(exons_dfs, psites_dfs)

    msg = "Finding all P-site intersections"
    logger.info(msg)

    sum_profiles = parallel.apply_parallel_iter(
        exons_psites,
        args.num_cpus,
        get_all_p_site_intersections,
        num_orfs,
        max_orf_len,
        progress_bar=True,
        total=len(exons_dfs),
        backend="multiprocessing",
    )

    msg = "Combining the ORF profiles into one matrix"
    logger.info(msg)

    def f(x, y):
        return x + y

    sum_profiles = functools.reduce(f, sum_profiles)
    sum_profiles_lil = sum_profiles.tolil()

    msg = "Flipping the reverse strand profiles"
    logger.info(msg)

    m_reverse = orfs["strand"] == "-"
    reverse_orfs = orfs[m_reverse]

    for idx, reverse_orf in tqdm.tqdm(reverse_orfs.iterrows()):
        orf_num = reverse_orf["orf_num"]

        if sum_profiles[orf_num].sum() == 0:
            continue

        orf_len = reverse_orf["orf_len"]
        dense = utils.to_dense(sum_profiles, orf_num, length=orf_len)
        dense = dense[::-1]
        sum_profiles_lil[orf_num, :orf_len] = dense

    msg = "Writing the sparse matrix to disk"
    logger.info(msg)
    math_utils.write_sparse_matrix(args.out, sum_profiles_lil)


if __name__ == "__main__":
    main()
