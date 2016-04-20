#! /usr/bin/env python3

import matplotlib
matplotlib.use('agg')

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

default_title = "Periodicity analysis"
default_xlabel = "Postion (bp). Translation start is at red line; end is at blue line; -12 is at green line."
default_ylabel = "Read count (starting at bp x)"
default_font_size = 15
default_series_label = ""
default_lengths = []

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script visualizes a provided metagene profile.")
    parser.add_argument('metagene_profile', help="The metagene profile (csv) file")
    parser.add_argument('length', help="The length to use for the visualization", type=int)
    parser.add_argument('out', help="The (output) image file")
    
    parser.add_argument('--title', help="The title for the figure", default=default_title)
    parser.add_argument('--xlabel', help="The label for the x-axis", default=default_xlabel)
    parser.add_argument('--ylabel', help="The label for the y-axis", default=default_ylabel)
    
    parser.add_argument('--series-label', help="The label for the legend", default=default_series_label)
    parser.add_argument('--font-size', help="The font size for the title, axis labels, and "
        "xticks labels", type=int, default=default_font_size)

    args = parser.parse_args()

    metagene_profile = pd.read_csv(args.metagene_profile)

    m_length = metagene_profile['length'] == args.length
    metagene_profile = metagene_profile[m_length]

    mask_starts = metagene_profile['type'] == 'start'
    start_counts = metagene_profile.loc[mask_starts, 'count'].values
    start_positions = metagene_profile.loc[mask_starts, 'position'].values

    # hackish, but a previous version took the positions as command line arguments
    args.start_upstream = start_positions[0]
    args.start_downstream = start_positions[-1]

    mask_ends = metagene_profile['type'] == 'end'
    end_counts = metagene_profile.loc[mask_ends, 'count'].values
    end_positions = metagene_profile.loc[mask_ends, 'position'].values

    args.end_upstream = end_positions[0]
    args.end_downstream = end_positions[-1]

    middle_counts = np.full(20, np.nan) # also some dummy counts to take space
    all_counts = np.concatenate([start_counts, middle_counts, end_counts])

    counts_sum = sum(start_counts) + sum(end_counts)
    scaled_counts = all_counts / counts_sum

    step = 10

    ymax = max(all_counts)
    scaled_ymax = max(scaled_counts)

    # first, get the labels for the bp around the translation start
    # this needs to be something like: range(-50, 21)
    # the +1 is to ensure we actually see the last start_downstream position in the labels
    xticklabels_start = range(args.start_upstream, args.start_downstream+1, step)
    xticks_start = range(0, len(start_counts), step)

    # now, the labels for the bp around the translation end site
    # this is also something like: range(-50, 21)
    xticklabels_end = range(args.end_upstream, args.end_downstream+1, step)
    xticks_end = range(len(start_counts) + len(middle_counts), len(all_counts), step)

    # join everything together
    xticklabels = np.concatenate([xticklabels_start, xticklabels_end])
    xticks = np.concatenate([xticks_start, xticks_end])

    fig, ax = plt.subplots(figsize=(10,5))

    # first, plot the counts
    ax.plot(all_counts, marker='o', color='g', label=args.series_label)

    # now, plot the lines to show translation start and end

    # start
    translation_start_position = -1*args.start_upstream
    ax.plot([translation_start_position, translation_start_position], [0,ymax], color='r', linestyle="--")

    # end
    translation_end_position = len(start_counts) + len(middle_counts) + -1*args.end_upstream
    ax.plot([translation_end_position, translation_end_position], [0,ymax], color='b', linestyle="--")

    # also show a line at -12
    minus_12 = translation_start_position - 12
    ax.plot([minus_12, minus_12], [0,ymax], color='g', linestyle="--")


    # next, draw some dots in the middle to highlight the gap
    # so, this should result in something like: [75, 80, 85]
    middle_dot_xvalues = range(len(start_counts) + int(step/2), len(start_counts) + len(middle_counts), int(step/2))
    middle_dot_yvalues = [ymax / 2] * len(middle_dot_xvalues)
    ax.scatter(middle_dot_xvalues, middle_dot_yvalues, color='k')

    # finally, add the labels, etc.
    ax.set_title(args.title, fontsize=args.font_size)
    ax.set_xlabel(args.xlabel, fontsize=args.font_size)
    ax.set_ylabel(args.ylabel, fontsize=args.font_size)
    ax.legend(loc="upper right")

    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, fontsize=args.font_size)

    ax.set_ylim((0, ymax))
    ax.set_xlim((0, len(all_counts)))

    fig.tight_layout()
    fig.savefig(args.out)

if __name__ == '__main__':
    main()

