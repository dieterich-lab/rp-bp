#! /usr/bin/env python3

import matplotlib
matplotlib.use('agg')

import argparse
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import pylab

import misc.utils as utils
import misc.stacked_bar_graph as stacked_bar_graph

# counting the transcriptome mapped reads does not work correctly right now...
#default_fields = ['raw_data_count', 'without_adapters_count', 'without_rrna_count', 'genome_count', 'unique_count', 'length_count', 'transcriptome_count']
#default_legend = ['poor_quality', 'rrna', 'no_genome_alignment', 'multimappers', 'wrong_length', 'no_transcript_alignment', 'usable_reads']

#without_rrna_fields = ['without_rrna_count', 'genome_count', 'unique_count', 'length_count', 'transcriptome_count']
#without_rrna_legend = ['no_genome_alignment', 'multimappers', 'wrong_length', 'no_transcript_alignment', 'usable_reads']

default_fields = ['raw_data_count', 'without_adapters_count', 'without_rrna_count', 'genome_count', 'unique_count', 'length_count']
default_legend = ['poor_quality', 'rrna', 'no_genome_alignment', 'multimappers', 'wrong_length', 'usable_reads']

without_rrna_fields = ['without_rrna_count', 'genome_count', 'unique_count', 'length_count']
without_rrna_legend = ['no_genome_alignment', 'multimappers', 'wrong_length', 'usable_reads']

default_title = "Filtered Reads"

default_font_size = 20
#default_ymax = 51000000
default_ystep = 10000000
default_ymax = 0
default_yticks = 6

# these come from the Spectral colormap
colors = [
    (0.61960786581039429, 0.0039215688593685627, 0.25882354378700256, 1.0),
    (0.91395617583218747, 0.3623990802203908, 0.27935410773052888, 1.0),
    (0.99346405267715454, 0.74771243333816528, 0.43529413143793744, 1.0),
    (0.99807766255210428, 0.99923106502084169, 0.74602077638401709, 1.0),
    (0.74771243333816528, 0.89803922176361084, 0.62745100259780884, 1.0),
    (0.3280276866520152, 0.68050751615973082, 0.68027683566598329, 1.0)
]


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script visualizes the number of reads lost at step in filtering. "
        "It contains several hard-coded values.")
    parser.add_argument('mapping_info', help="The (csv) mapping information file")
    parser.add_argument('out', help="The output image file")
    parser.add_argument('--title', help="The title for the plot", default=default_title)
    parser.add_argument('--fields', help="The fields to use from the mapping info file", 
        nargs='+', default=default_fields)
    parser.add_argument('--legend', help="The text to use for the kept fields in the "
        "images. The order should match the order of fields.", nargs='+', default=default_legend)
    parser.add_argument('--without-rrna', help="If this flag is present, then a default set "
        "of fields excluding the reads mapping to ribosomal sequences will be used.", 
        action='store_true')
    parser.add_argument('--font-size', help="The font size to use for everything in "
        "the plot", type=int, default=default_font_size)
    parser.add_argument('--ymax', help="If present, this value will be used as the maximum for "
        "the y-axis.", type=int, default=default_ymax)
    parser.add_argument('--ystep', help="If --ymax is specified, then this value will be used as "
        "the step size for labels on the y-axis.", type=int, default=default_ystep)

    parser.add_argument('--default-yticks', help="If --ymax and --ystep are not specified, "
        "then this many ticks will be used", type=int, default=default_yticks)
    
    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    if args.without_rrna:
        args.legend = without_rrna_legend
        args.fields = without_rrna_fields

    # flip the legend so that the labels and fields match
    args.legend = args.legend[::-1]

    # read the data
    mapping_info = pd.read_csv(args.mapping_info, comment='#')
    mapping_info_viz = mapping_info[args.fields]
    data_np = np.array(mapping_info_viz)

    # add an extra column so the diff counts will work
    zeros = np.zeros((data_np.shape[0], 1))
    data_np = np.append(data_np, zeros, axis=1)

    # get the diffs so the stacks work correctly
    # this is probably somehow redundante, but necessary for stackedBarGraph
    diff = np.array([data_np[:,i] - data_np[:,i+1] for i in range(data_np.shape[1]-2, -1, -1)])
    diff = diff.T

    logging.debug(data_np)
    logging.debug(diff)

    # we can use the note column to get the sample names
    names = mapping_info['note']

    # use the stackedBarPlot class to plot
    SBG = stacked_bar_graph.StackedBarGrapher()
    fig, ax = plt.subplots()
    gap = 0

    if args.ymax != default_ymax:
        ylabels = np.arange(0, args.ymax, args.ystep)
        yticks = (ylabels, ylabels)
    else:
        yticks = args.default_yticks

    #edge_colors = [matplotlib.colors.cnames['darkgrey']] * (len(names)+1)
    edge_colors = None

    bars = SBG.stackedBarPlot(ax,
                        diff,
                        colors,
                        edgeCols=edge_colors,
                        xLabels=names,
                        yTicks = yticks,
                        gap=gap,
                        endGaps=True,
                        legend=args.legend,
                        xlabel='Samples',
                        ylabel='Reads',
                        log=False,
                        legend_loc='upper center',
                        legend_bbox_to_anchor=(0.5, -0.6),
                        legend_ncol=2,
                        font_size=args.font_size
                        )
    lgd = ax.get_legend()
    ltext = lgd.get_texts()
    plt.setp(ltext, fontsize=20)
    lgd.set_title("Filter reason")
    plt.setp(lgd.get_title(), fontsize=args.font_size)
    
    ax.set_title(args.title)
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

    ax.tick_params(axis='both', which='major', labelsize=args.font_size)

    fig.tight_layout()
    fig.savefig(args.out, bbox_inches='tight', bbox_extra_artist=[lgd])
    #fig.savefig(args.out, additional_artists = art)

if __name__ == '__main__':
    main()

