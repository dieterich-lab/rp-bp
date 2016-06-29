#! /usr/bin/env python3

import matplotlib
matplotlib.use('agg')

import argparse
import logging

import matplotlib.pyplot as plt
import numpy as np
import scipy.io

import misc.bio as bio
import misc.parallel as parallel
import misc.utils as utils

import riboutils.ribo_filenames as filenames

default_image_type = 'eps'
default_num_cpus = 1

def get_windows(profile):
    
    profile = profile / np.max(profile)
    
    orf_len = len(profile)
    if orf_len < 42:
        # we would return first window and exit
        first_window = profile[:21]
        return (first_window, None, None)

    first_window, middle_window, last_window = np.split(profile, [21, orf_len-21])

    # now, pull together and sum up all intermediate windows (of length 21)
    # cheat a bit, and just split split the middle into 21-bp windows, drop the last window
    indices = np.arange(21, len(middle_window), 21)
    middle_windows = np.split(middle_window, indices)[:-1]
    
    return first_window, middle_windows, last_window

def get_profile(orf, profiles):
    orf_num = orf['orf_num']
    orf_len = orf['orf_len']

    if orf_len < 21:
        return None

    profile = utils.to_dense(profiles, orf_num, length=orf_len)
    return profile

def plot_windows(windows, title, out):

    windows_np = np.array(windows)
    first_windows = windows_np[:,0]

    last_windows = windows_np[:,2] 
    last_windows = np.array([lw for lw in last_windows if lw is not None])

    middle_windows = windows_np[:,1] 
    middle_windows = [mw for mw in middle_windows if mw is not None]
    middle_windows = utils.flatten_lists(middle_windows)
    middle_windows = np.array(middle_windows)

    ind = np.arange(21)  # the x locations for the groups
    width = 0.5       # the width of the bars

    fig, axes = plt.subplots(ncols=3, sharey=True, sharex=True, figsize=(10,5))

    # the first window
    first_means = np.mean(first_windows, axis=0)
    first_var = np.var(first_windows, axis=0)
    rects_first = axes[0].bar(ind, first_means, width, color='g', yerr=first_var)

    # the middle windows
    middle_means = np.mean(middle_windows, axis=0)
    middle_var = np.var(middle_windows, axis=0)
    rects_middle = axes[1].bar(ind, middle_means, width, color='g', yerr=middle_var)

    # the last window
    last_means = np.mean(last_windows, axis=0)
    last_var = np.var(last_windows, axis=0)
    rects_last = axes[2].bar(ind, last_means, width, color='g', yerr=last_var)

    axes[0].set_xlim((-width, 21))
    axes[0].set_ylim((0, 1.05))
   
    axes[0].set_title('First 21-bp window')
    axes[1].set_title('All 21-bp windows in middle')
    axes[2].set_title('Last 21-bp window')
    fig.suptitle(title)

    fig.savefig(out, bbox_inches='tight')

def extract_profiles_and_plot(g, profiles, args):
    orf_type = g['orf_type'].iloc[0]

    msg = "ORF type: {}".format(orf_type)
    logging.info(msg)

    msg = "Extracting profiles"
    logging.debug(msg)
    g_profiles = parallel.apply_df_simple(g, get_profile, profiles, progress_bar=True)

    msg = "Slicing the profiles into windows"
    logging.debug(msg)
    windows = parallel.apply_parallel_iter(g_profiles, args.num_cpus, get_windows, progress_bar=True)
    
    msg = "Plotting the profile statistics"
    logging.debug(msg)

    out = filenames.get_orf_type_profile_image(args.out, orf_type, args.image_type)
    plot_windows(windows, orf_type, out)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script visualizes the metagene profiles for each ORF type "
        "present in a given BED12+ file. It visualizes the mean and variance of normalized "
        "profiles in the first 21-bp, last 21-bp, and across all other 21-bp windows.")

    parser.add_argument('orfs', help="The BED12+ file containing the ORFs")
    parser.add_argument('profiles', help="The (mtx) file containing the ORF profiles")
    parser.add_argument('out', help="The base output name. The output filenames will be of "
        "the form: <out>.<orf-type>.<image-type>.")

    parser.add_argument('--image-type', help="The type of image files to create. The type "
        "must be recognized by matplotlib.", default=default_image_type)

    parser.add_argument('--num-cpus', help="The number of cores to use", type=int,
        default=default_num_cpus)
    
    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    msg = "Reading ORFs"
    logging.info(msg)
    orfs = bio.read_bed(args.orfs)

    msg = "Reading profiles"
    logging.info(msg)
    profiles = scipy.io.mmread(args.profiles).tocsr()

    msg = "Extracting the metagene profiles and creating the images"
    logging.info(msg)

    orf_type_groups = orfs.groupby('orf_type')
    orf_type_groups.apply(extract_profiles_and_plot, profiles, args)

    msg = "Finished"
    logging.info(msg)

if __name__ == '__main__':
    main()
