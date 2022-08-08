#! /usr/bin/env python3

import matplotlib

matplotlib.use("agg")

import argparse
import yaml
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import seaborn as sns

sns.set(style="white")

import pbio.misc.mpl_utils as mpl_utils
import pbio.misc.logging_utils as logging_utils

import pbio.ribo.ribo_utils as ribo_utils


logger = logging.getLogger(__name__)

default_fontsize = 20
default_legend_fontsize = 15

default_ymax = None  # 1e8+1
default_ystep = 1e7

default_alignment_counts_order = [
    "raw_data_count",
    "without_adapters_count",
    "without_rrna_count",
    "genome_count",
    "unique_count",
    "length_count",
]

default_alignment_counts_names = [
    "Poor quality",
    "Ribosomal",
    "No alignment",
    "Multimappers",
    "Non-periodic",
    "Usable",
]

without_rrna_order = [
    "without_rrna_count",
    "genome_count",
    "unique_count",
    "length_count",
]

without_rrna_names = ["No alignment", "Multimappers", "Non-periodic", "Usable"]


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Visualize the read counts at each filtering step.",
    )

    parser.add_argument(
        "alignment_counts",
        help="The (csv) alignment counts " "(created with get-all-filtering-counts)",
    )
    parser.add_argument("out", help="The output image file")

    parser.add_argument(
        "--alignment-counts-order",
        help="The fields to use "
        "from the alignment_counts file. The order should range from the "
        "least strict filter to most strict filter.",
        nargs="+",
        default=default_alignment_counts_order,
    )

    parser.add_argument(
        "--alignment-counts-names",
        help="The text to use for "
        "the kept fields in the plot. The order of the names must match that "
        "of --alignment-counts-order.",
        nargs="+",
        default=default_alignment_counts_names,
    )

    parser.add_argument(
        "--without-rrna",
        help="If this flag is present, then "
        "a default set of fields excluding the reads mapping to ribosomal "
        "sequences will be used.",
        action="store_true",
    )

    parser.add_argument(
        "--config",
        help="""The config file, if using
        pretty names, "riboseq_sample_name_map" must be defined""",
        type=str,
        default=None,
    )

    parser.add_argument("--title", help="The title of the plot", default=None)

    parser.add_argument(
        "--fontsize",
        help="The font size to use for most of " "the text in the plot",
        type=int,
        default=default_fontsize,
    )

    parser.add_argument(
        "--legend-fontsize",
        help="The font size to use for " "the legend labels",
        type=int,
        default=default_legend_fontsize,
    )

    parser.add_argument(
        "--ymax", help="The maximum for the y-axis", type=int, default=default_ymax
    )
    parser.add_argument(
        "--ystep",
        help="The step size for ticks on the " "y-axis.",
        type=int,
        default=default_ystep,
    )

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    if args.without_rrna:
        msg = "Using the default without rrna field order"
        logger.info(msg)

        args.alignment_counts_order = without_rrna_order
        args.alignment_counts_names = without_rrna_names

    # we really need the fields from most- to least-restrictive
    args.alignment_counts_order = args.alignment_counts_order[::-1]
    args.alignment_counts_names = args.alignment_counts_names[::-1]

    msg = "Reading counts"
    logger.info(msg)

    alignment_counts = pd.read_csv(args.alignment_counts)
    alignment_counts = alignment_counts.sort_values("note")

    msg = "Calculating the diff counts"
    logger.info(msg)
    alignment_diff_counts = mpl_utils.get_diff_counts(
        alignment_counts[args.alignment_counts_order]
    )

    df = pd.DataFrame(alignment_diff_counts)
    df.columns = args.alignment_counts_names

    names = alignment_counts["note"].reset_index(drop=True)
    df["name"] = names

    if args.config:
        try:
            config = yaml.load(open(args.config), Loader=yaml.FullLoader)
            sample_name_map = ribo_utils.get_sample_name_map(config)
            df["display_name"] = df["name"].apply(lambda x: sample_name_map[x])
        except:
            msg = 'Fall back to "name", cannot fetch "display_name" from config file.'
            logger.warning(msg)
            df["display_name"] = df["name"]
    else:
        df["display_name"] = df["name"]

    msg = "Creating the stacked bar chart"
    logger.info(msg)

    fig, ax = plt.subplots()

    pal = sns.palettes.color_palette(
        palette="Set3", n_colors=len(args.alignment_counts_names)
    )

    gap = 0.15

    # if we aren't given information about the y-axis, try to guess
    if args.ymax is None:
        field = args.alignment_counts_order[-2]
        max_count = alignment_counts[field].max()

        if args.ystep > max_count:
            args.ystep = np.ceil(max_count / 4)

        args.ymax = (np.ceil(max_count / args.ystep) * args.ystep) + 1

    yticks = np.arange(0, args.ymax, args.ystep)

    bars = mpl_utils.create_stacked_bar_graph(
        ax,
        alignment_diff_counts,
        colors=pal,
        x_tick_labels=df["display_name"],
        y_ticks=yticks,
        y_tick_labels=yticks,
        gap=gap,
        end_gaps=True,
        stack_labels=args.alignment_counts_names,
        y_title="Reads",
        log=False,
        font_size=args.fontsize,
        edge_colors="0.5",
    )

    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.6),
        ncol=3,
        fontsize=args.legend_fontsize,
        title="Filter",
        frameon=True,
    )

    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter("%.0e"))
    mpl_utils.set_label_fontsize(ax, args.fontsize)
    mpl_utils.set_legend_title_fontsize(ax, args.fontsize)

    if args.title is not None:
        ax.set_title(args.title, fontsize=args.fontsize)

    msg = "Writing the plot to disk"
    logger.info(msg)
    fig.savefig(args.out, bbox_inches="tight")


if __name__ == "__main__":
    main()
