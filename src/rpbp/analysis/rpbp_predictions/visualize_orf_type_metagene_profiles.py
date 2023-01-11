#! /usr/bin/env python3

import matplotlib

import argparse
import logging
import sys
import yaml

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io
import tqdm

import pbiotools.utils.bed_utils as bed_utils
import pbiotools.misc.parallel as parallel
import pbiotools.misc.logging_utils as logging_utils
import pbiotools.misc.utils as utils

import rpbp.ribo_utils.filenames as filenames

matplotlib.use("agg")

logger = logging.getLogger(__name__)

default_image_type = "eps"
default_title = ""
default_min_profile = 5
default_max_orfs = 10000


def get_windows(profile):

    profile = profile / np.max(profile)

    orf_len = len(profile)
    if orf_len < 42:
        # we would return first window and exit
        first_window = profile[:21]
        return (first_window, None, None)

    first_window, middle_window, last_window = np.split(profile, [21, orf_len - 21])

    # now, pull together and sum up all intermediate windows (of length 21)
    # cheat a bit, and just split split the middle into 21-bp windows, drop the last window
    indices = np.arange(21, len(middle_window), 21)
    middle_windows = np.split(middle_window, indices)[:-1]

    return first_window, middle_windows, last_window


def get_profile(orf, profiles, min_profile):
    orf_num = orf["orf_num"]
    orf_len = orf["orf_len"]

    if orf_len < 21:
        return None

    profile = utils.to_dense(profiles, orf_num, length=orf_len)

    if sum(profile) < min_profile:
        return None

    return profile


def plot_windows(windows, title, out):

    if len(windows) == 0:
        msg = "Did not find any windows for: {}".format(title)
        logger.warning(msg)
        return

    windows_np = np.array(windows)
    first_windows = windows_np[:, 0]

    last_windows = windows_np[:, 2]
    last_windows = np.array([lw for lw in last_windows if lw is not None])

    middle_windows = windows_np[:, 1]
    middle_windows = [mw for mw in middle_windows if mw is not None]
    middle_windows = utils.flatten_lists(middle_windows)
    middle_windows = np.array(middle_windows)

    if (len(last_windows) == 0) or (len(middle_windows) == 0):
        msg = "No long ORFs for: {}".format(title)
        logger.warning(msg)
        return

    ind = np.arange(21)  # the x locations for the groups
    width = 0.5  # the width of the bars

    fig, axes = plt.subplots(ncols=3, sharey=True, sharex=True, figsize=(10, 5))

    # the first window
    first_means = np.mean(first_windows, axis=0)
    first_var = np.var(first_windows, axis=0)
    rects_first = axes[0].bar(ind, first_means, width, color="g", yerr=first_var)

    # the middle windows
    middle_means = np.mean(middle_windows, axis=0)
    middle_var = np.var(middle_windows, axis=0)
    rects_middle = axes[1].bar(ind, middle_means, width, color="g", yerr=middle_var)

    # the last window
    last_means = np.mean(last_windows, axis=0)
    last_var = np.var(last_windows, axis=0)
    rects_last = axes[2].bar(ind, last_means, width, color="g", yerr=last_var)

    axes[0].set_xlim((-width, 21))
    axes[0].set_ylim((0, 1.05))

    axes[0].set_title("First 21-bp window")
    axes[1].set_title("All 21-bp windows in middle")
    axes[2].set_title("Last 21-bp window")
    fig.suptitle(title)

    msg = "Saving figure to: {}".format(out)
    logger.debug(msg)
    fig.savefig(out, bbox_inches="tight")


def extract_profiles_and_plot_strand(g, profiles, orf_type, strand, args):
    m_strand = g["strand"] == strand

    msg = "Extracting profiles"
    logger.debug(msg)

    # we will manually build up the list so we can quit as soon as possible
    df = g[m_strand]

    # quickly figure out which profiles have a sum we can use
    orf_nums = np.array(df["orf_num"])

    # we may have only calculated profiles for a few ORFs (using --num-orfs)
    # so do not keep indices which are too large
    num_profiles = profiles.shape[0]
    m_orfs = orf_nums < num_profiles
    orf_nums = orf_nums[m_orfs]

    g_profiles = profiles[orf_nums]
    s = g_profiles.sum(axis=1)
    m_sum = np.array(s > args.min_profile)[:, 0]

    # we need to extend m_sum so it is the correct length
    num_orfs = len(orf_nums)
    diff = len(m_orfs) - num_orfs
    m_orfs = np.concatenate([m_sum, np.array([False] * diff)]).astype(bool)

    g_profiles = []
    df_rows = tqdm.tqdm(
        df[m_orfs].iterrows(), total=len(df[m_orfs]), leave=True, file=sys.stdout
    )

    for row in df_rows:
        profile = get_profile(row[1], profiles, args.min_profile)
        if profile is not None:
            g_profiles.append(profile)

            if len(g_profiles) >= args.max_orfs:
                break

    msg = "Slicing the profiles into windows"
    logger.debug(msg)
    windows = parallel.apply_iter_simple(g_profiles, get_windows, progress_bar=True)

    msg = "Plotting the profile statistics"
    logger.debug(msg)

    out = filenames.get_orf_type_profile_image(
        args.out, orf_type, strand, args.image_type
    )

    title = "{}: {}, strand: {} ({})".format(
        args.title, orf_type, strand, len(df[m_orfs])
    )
    plot_windows(windows, title, out)


def extract_profiles_and_plot(g, profiles, args):
    orf_type = g["orf_type"].iloc[0]

    msg = "ORF type: {}".format(orf_type)
    logger.info(msg)

    msg = "Strand: -"
    logger.info(msg)

    strand = "-"
    extract_profiles_and_plot_strand(g, profiles, orf_type, strand, args)

    msg = "Strand: +"
    logger.info(msg)

    strand = "+"
    extract_profiles_and_plot_strand(g, profiles, orf_type, strand, args)


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script visualizes the metagene profiles for each ORF type "
        "present in a given BED12+ file. It visualizes the mean and variance of normalized "
        "profiles in the first 21-bp, last 21-bp, and across all other 21-bp windows.",
    )

    parser.add_argument("config", help="The (yaml) config file")
    parser.add_argument("orfs", help="The BED12+ file containing the ORFs")
    parser.add_argument("profiles", help="The (mtx) file containing the ORF profiles")
    parser.add_argument(
        "out",
        help="The base output name. The output filenames will be of "
        "the form: <out>.<orf-type>.<image-type>.",
    )

    parser.add_argument(
        "--min-profile",
        help="The minimum value of the sum over the profile "
        "to include it in the analysis",
        type=float,
        default=default_min_profile,
    )

    parser.add_argument(
        "--max-orfs",
        help="At most this many ORFs of each type will be "
        "used to create the figures. They will be sampled randomly from among those "
        "which meet the min-profile constraint.",
        type=int,
        default=default_max_orfs,
    )

    parser.add_argument(
        "--title",
        help="The prefix to use for the title of the plots",
        default=default_title,
    )

    parser.add_argument(
        "--image-type",
        help="The type of image files to create. The type "
        "must be recognized by matplotlib.",
        default=default_image_type,
    )

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    config = yaml.load(open(args.config), Loader=yaml.FullLoader)
    note = config.get("note", None)

    msg = "Reading ORFs"
    logger.info(msg)
    orfs = bed_utils.read_bed(args.orfs)

    msg = "Adding labels to ORFs"
    logger.info(msg)

    labeled_orfs = filenames.get_labels(
        config["genome_base_path"], config["genome_name"], note=note
    )
    cols = ["id", "orf_type"]
    labels_df = bed_utils.read_bed(labeled_orfs)[cols]
    orfs = pd.merge(orfs, labels_df, on="id", how="left")

    msg = "Reading profiles"
    logger.info(msg)
    profiles = scipy.io.mmread(args.profiles).tocsr()

    msg = "Extracting the metagene profiles and creating the images"
    logger.info(msg)

    orf_type_groups = orfs.groupby("orf_type")
    orf_type_groups.apply(extract_profiles_and_plot, profiles, args)

    msg = "Finished"
    logger.info(msg)


if __name__ == "__main__":
    main()
