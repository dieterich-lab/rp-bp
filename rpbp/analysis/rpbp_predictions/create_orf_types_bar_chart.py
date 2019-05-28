#! /usr/bin/env python3

import matplotlib
matplotlib.use('agg')

import argparse

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns; sns.set(style='white')

import pbio.utils.bed_utils as bed_utils
import pbio.misc.mpl_utils as mpl_utils
import pbio.ribo.ribo_utils as ribo_utils

import logging
import pbio.misc.logging_utils as logging_utils
logger = logging.getLogger(__name__)

default_title = None
default_fontsize = 20
default_legend_fontsize = 15
default_ymax = 1e4

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script creates a bar chart which shows the count of "
        "each ORF type in a given BED12+ file. Optionally, the ORFs can be "
        "grouped into similar types.")

    parser.add_argument('orfs', help="The BED12+ file with the ORFs")
    parser.add_argument('out', help="The output (image) file")

    parser.add_argument('--title', help="The title to use for the plot", 
        default=default_title)

    parser.add_argument('--use-groups', help="If this flag is given, the ORFs "
        "will be grouped", action='store_true')

    parser.add_argument('--fontsize', default=default_fontsize)
    parser.add_argument('--legend-fontsize', default=default_legend_fontsize)
    parser.add_argument('--ymax', type=int, default=default_ymax)
    
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()

    msg = "Reading bed file"
    logger.info(msg)

    bed = bed_utils.read_bed(args.orfs)

    if args.use_groups:
        bed['orf_type_group'] = bed['orf_type'].map(
            ribo_utils.orf_type_labels_reverse_mapping)

        orf_type_counts = bed.groupby(['orf_type_group', 'strand']).size()
        orf_type_counts = orf_type_counts.reset_index(name="count")
        orf_type_counts['display_name'] = orf_type_counts['orf_type_group'].map(
            ribo_utils.orf_type_labels_display_name_map)
    else:
        orf_type_counts = bed.groupby(['orf_type', 'strand']).size()
        orf_type_counts = orf_type_counts.reset_index(name="count")
        orf_type_counts['display_name'] = orf_type_counts['orf_type'].map(
            ribo_utils.orf_type_display_name_map)

    msg = "Creating the bar chart"

    color = sns.palettes.color_palette("Set3", n_colors=3)

    fig, ax = plt.subplots(figsize=(9,5))
    sns.barplot(
        x="display_name",
        y="count",
        hue="strand",
        data=orf_type_counts,
        ax=ax,
        zorder=-1,
        palette='Set3',
        log=True
    )

    sns.despine()

    ax.legend(
        loc='upper right', 
        bbox_to_anchor=(1.0, 0.95), 
        fontsize=args.legend_fontsize, 
        frameon=True, 
        framealpha=0.9,
        title="Strand"
    )
    mpl_utils.set_legend_title_fontsize(ax, args.fontsize)

    #ax.set_yscale('log')
    #ax.set_ylim((1, args.ymax))

    ax.set_ylabel("Number of ORFs", fontsize=args.fontsize)
    ax.set_xlabel("", fontsize=0)

    # rotate the ORF type names
    mpl_utils.set_ticklabels_fontsize(ax, args.fontsize)
    mpl_utils.set_ticklabel_rotation(ax, axis='x', rotation=90)

    # place the ORF type names in the middle of the bar
    for ticklabel in ax.xaxis.get_ticklabels():    
        p = ticklabel.get_position()
        ticklabel.set_position((p[0], 0.1))
        ticklabel.set_verticalalignment('bottom')
        
    if args.title is not None:
        ax.set_title(args.title, fontsize=args.fontsize)
        
    if args.out is not None:
        fig.savefig(args.out, bbox_inches='tight')

if __name__ == '__main__':
    main()
