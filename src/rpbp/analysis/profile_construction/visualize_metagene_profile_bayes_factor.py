#! /usr/bin/env python3

import matplotlib

import argparse
import matplotlib.pyplot as plt
import pandas as pd

matplotlib.use("agg")

default_title = "Metagene profile Bayes' factors"
default_xlabel = "Offset, relative to translation \ninitiation site"
default_ylabel = "Bayes' factor"
default_font_size = 15

default_series_label = ""


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script visualizes the Bayes' factors for a metagene profile.\n\n"
        "This script contains some hard-coded field names.",
    )
    parser.add_argument("bayes_factors", help="The metagene profile (csv) file")
    parser.add_argument("length", help="The profile lengths to visualize", type=int)
    parser.add_argument("out", help="The (output) image file")

    parser.add_argument(
        "--title", help="The title for the figure", default=default_title
    )
    parser.add_argument(
        "--xlabel", help="The label for the x-axis", default=default_xlabel
    )
    parser.add_argument(
        "--ylabel", help="The label for the y-axis", default=default_ylabel
    )
    parser.add_argument(
        "--series-label", help="The label for the legend", default=default_series_label
    )
    parser.add_argument(
        "--font-size",
        help="The font size for the title, axis labels, and " "xticks labels",
        type=int,
        default=default_font_size,
    )

    args = parser.parse_args()

    bayes_factors = pd.read_csv(args.bayes_factors)

    mask_length = bayes_factors["length"] == args.length
    group = bayes_factors.loc[mask_length]

    bfs = group["bayes_factor_mean"]
    offsets = group["offset"]
    bf_range = max(bfs) - min(bfs)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(offsets, bfs, label=args.series_label, color="b")
    ax.scatter(offsets, bfs, color="b")

    xlim = (min(offsets), max(offsets))

    ymin = min(bfs) - 0.1 * bf_range
    ymax = max(bfs) + 0.1 * bf_range
    ylim = (ymin, ymax)

    # and draw a line at "bf=5"
    plt.plot(xlim, (5, 5), color="k", linewidth=2, linestyle=":")

    # and a horizontal line at the maximum bf
    plt.plot(xlim, (max(bfs), max(bfs)), color="r", linewidth=1, linestyle="-.")

    # and a vertical line at "offset=-12"
    ax.plot((-12, -12), ylim, color="g", linestyle="--")

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    # finally, add the labels, etc.
    plt.suptitle(args.title, fontsize=args.font_size, y=1.03)
    ax.set_xlabel(args.xlabel, fontsize=args.font_size)
    ax.set_ylabel(args.ylabel, fontsize=args.font_size)

    ax.tick_params(axis="both", which="major", labelsize=args.font_size)
    # ax.legend(loc="upper right")

    fig.tight_layout()
    fig.savefig(args.out, bbox_inches="tight")


if __name__ == "__main__":
    main()
