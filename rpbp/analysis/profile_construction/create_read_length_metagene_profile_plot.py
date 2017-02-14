#! /usr/bin/env python3

import matplotlib
matplotlib.use('agg')

import argparse
import matplotlib.pyplot as plt
import seaborn as sns; sns.set(style='white')
import numpy as np
import pandas as pd

default_title = "Periodicity analysis"
default_xlabel_start = "Position of P-site relative to start (nt)\nRed: TIS. Green: TIS -12"
default_xlabel_end = "Position of P-site relative to stop (nt)\nBlue: Translation termination"
default_ylabel = "Read count (starting at bp x)"
default_font_size = 15
default_series_label = ""
default_lengths = []

default_start_upstream = -50
default_start_downstream = 21

default_end_upstream = -50
default_end_downstream = 21

default_step = 10

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script visualizes a provided metagene profile.")
    parser.add_argument('metagene_profile', help="The metagene profile (csv) file")
    parser.add_argument('length', help="The length to use for the visualization", type=int)
    parser.add_argument('out', help="The (output) image file")
    
    parser.add_argument('--title', help="The title for the figure", default=default_title)
    parser.add_argument('--xlabel-start', help="The label for the x-axis for TIS", default=default_xlabel_start)
    parser.add_argument('--xlabel-end', help="The label for the x-axis for TTS", default=default_xlabel_end)
    parser.add_argument('--ylabel', help="The label for the y-axis", default=default_ylabel)

    parser.add_argument('--step', help="The step for the x-axis", type=int, default=default_step)
    
    parser.add_argument('--series-label', help="The label for the legend", default=default_series_label)
    parser.add_argument('--font-size', help="The font size for the title, axis labels, and "
        "xticks labels", type=int, default=default_font_size)

    parser.add_argument('--start-upstream', type=int, help="The number of bp "
        "upstream of the TIS to show. This value should be negative.", 
        default=default_start_upstream)
    parser.add_argument('--start-downstream', type=int, help="The number of "
        "bp downstream of the TIS to show", default=default_start_downstream)
        
    parser.add_argument('--end-upstream', type=int, help="The number of bp "
        "upstream of the TTS to show. This value should be negative.", 
        default=default_end_upstream)
    parser.add_argument('--end-downstream', type=int, help="The number of "
        "bp downstream of the TTS to show", default=default_end_downstream)

    parser.add_argument('--use-entire-profile', help="If this flag is given, "
        "the the start_upstream, etc., values will be taken as everything in "
        "provided profile rather than the command line parameters",
        action='store_true')

    args = parser.parse_args()

    metagene_profile = pd.read_csv(args.metagene_profile)

    m_length = metagene_profile['length'] == args.length
    metagene_profile = metagene_profile[m_length]

    mask_starts = metagene_profile['type'] == 'start'
    mask_ends = metagene_profile['type'] == 'end'
    
    if args.use_entire_profile:
        args.start_upstream = start_positions[0]
        args.start_downstream = start_positions[-1]
    
        args.end_upstream = end_positions[0]
        args.end_downstream = end_positions[-1]

    else:
        m_start_upstream = metagene_profile['position'] >= args.start_upstream
        m_start_downstream = metagene_profile['position'] <= args.start_downstream
        
        m_end_upstream = metagene_profile['position'] >= args.end_upstream
        m_end_downstream = metagene_profile['position'] <= args.end_downstream

        mask_starts = m_start_upstream & m_start_downstream & mask_starts
        mask_ends = m_end_upstream & m_end_downstream & mask_ends
    
    start_counts = metagene_profile.loc[mask_starts, 'count'].values
    start_positions = metagene_profile.loc[mask_starts, 'position'].values

    end_counts = metagene_profile.loc[mask_ends, 'count'].values
    end_positions = metagene_profile.loc[mask_ends, 'position'].values

    ymax = max(max(start_counts), max(end_counts))
    ymax = 1.1 * ymax

    # first, get the labels for the bp around the translation start
    # this needs to be something like: range(-50, 21)
    # the +1 is to ensure we actually see the last start_downstream position in the labels
    xticklabels_start = range(args.start_upstream, args.start_downstream+1, args.step)
    xticks_start = range(0, len(start_counts), args.step)

    # now, the labels for the bp around the translation end site
    # this is also something like: range(-50, 21)
    xticklabels_end = range(args.end_upstream, args.end_downstream+1, args.step)
    #xticks_end = range(len(start_counts) + len(middle_counts), len(all_counts), args.step)
    xticks_end = range(0, len(end_counts), args.step)

    fig, axes = plt.subplots(ncols=2, figsize=(10,5), sharey=True)
    colors = sns.color_palette("Set3", 3)

    # first, the start counts
    x_1 = start_counts[0::3]
    x_2 = start_counts[1::3]
    x_3 = start_counts[2::3]

    x_1_pos = np.arange(len(start_counts))[0::3]
    x_2_pos = np.arange(len(start_counts))[1::3]
    x_3_pos = np.arange(len(start_counts))[2::3]

    x_1_rects = axes[0].bar(x_1_pos, x_1, width=1, color=colors[0])
    x_2_rects = axes[0].bar(x_2_pos, x_2, width=1, color=colors[1])
    x_3_rects = axes[0].bar(x_3_pos, x_3, width=1, color=colors[2])

    # now, the end counts
    
    x_1 = end_counts[0::3]
    x_2 = end_counts[1::3]
    x_3 = end_counts[2::3]

    x_1_pos = np.arange(len(end_counts))[0::3]
    x_2_pos = np.arange(len(end_counts))[1::3]
    x_3_pos = np.arange(len(end_counts))[2::3]

    x_1_rects = axes[1].bar(x_1_pos, x_1, width=1, color=colors[0])
    x_2_rects = axes[1].bar(x_2_pos, x_2, width=1, color=colors[1])
    x_3_rects = axes[1].bar(x_3_pos, x_3, width=1, color=colors[2])


    # now, plot the lines to show translation start and end

    # start
    translation_start_position = -1*args.start_upstream
    axes[0].plot([translation_start_position, translation_start_position], [0,ymax], color='r', linestyle="--")
    
    # also show a line at -12
    minus_12 = translation_start_position - 12
    axes[0].plot([minus_12, minus_12], [0,ymax], color='g', linestyle="--")

    # end
    translation_end_position = -1 * args.end_upstream
    axes[1].plot([translation_end_position, translation_end_position], [0,ymax], color='b', linestyle="--")

    # finally, add the labels, etc.
    axes[0].set_ylabel(args.ylabel, fontsize=args.font_size)

    axes[0].set_xlabel(args.xlabel_start, fontsize=args.font_size)
    axes[1].set_xlabel(args.xlabel_end, fontsize=args.font_size)

    if args.title is not None:
        fig.suptitle(args.title, fontsize=args.font_size, y=1.03)

    axes[0].set_xticks(xticks_start)
    axes[0].set_xticklabels(xticklabels_start, fontsize=args.font_size)
    
    axes[1].set_xticks(xticks_end)
    axes[1].set_xticklabels(xticklabels_end, fontsize=args.font_size)

    axes[0].set_ylim((0, ymax))

    # hide the spines between ax and ax2
    axes[0].spines['right'].set_visible(False)
    axes[1].spines['left'].set_visible(False)
    axes[0].yaxis.tick_left()
    axes[0].tick_params(labelright='off')
    axes[1].yaxis.tick_right()
    
    d = .0075 # how big to make the diagonal lines in axes coordinates
    vertical_adj = 3 # a factor for controlling the vertical slope
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=axes[0].transAxes, color='k', clip_on=False)
    axes[0].plot((1-d,1+d), (-vertical_adj*d,+vertical_adj*d), **kwargs)
    axes[0].plot((1-d,1+d),(1-vertical_adj*d,1+vertical_adj*d), **kwargs)

    kwargs.update(transform=axes[1].transAxes)  # switch to the bottom axes
    axes[1].plot((-d,+d), (1-vertical_adj*d,1+vertical_adj*d), **kwargs)
    axes[1].plot((-d,+d), (-vertical_adj*d,+vertical_adj*d), **kwargs)

    fig.tight_layout()
    fig.savefig(args.out, bbox_inches='tight')

if __name__ == '__main__':
    main()

