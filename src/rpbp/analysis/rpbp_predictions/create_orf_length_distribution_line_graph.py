#! /usr/bin/env python3

import matplotlib

matplotlib.use("agg")
matplotlib.rc("text", usetex=True)

import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.stats

import pbio.utils.bio as bio
import pbio.utils.bed_utils as bed_utils
import pbio.misc.latex as latex
import pbio.misc.math_utils as math_utils

import pbio.ribo.ribo_utils as ribo_utils

default_uniprot = ""
default_uniprot_label = "UniRef"
default_title = ""


def get_orf_lengths(orfs, orf_types):
    m_orf_type = orfs["orf_type"].isin(orf_types)
    lengths = np.array(orfs.loc[m_orf_type, "orf_len"])
    return lengths


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script creates a line graph showing the length distributions "
        "of the various types of ORFs. Optionally, it can also include the length "
        "distribution of ORFs downloaded from uniprot. If uniprot ORFs are given, then the "
        "KL-divergence between the type distributions and the uniprot ORFs is calculated.",
    )

    parser.add_argument("orfs", help="The BED12+ file with the ORFs")
    parser.add_argument("out", help="The output (image) file")

    parser.add_argument(
        "--uniprot",
        help="The uniprot ORF lengths, if available",
        default=default_uniprot,
    )
    parser.add_argument(
        "--uniprot-label",
        help="The label to use for the uniprot ORFs in " "the plot",
        default=default_uniprot_label,
    )
    parser.add_argument(
        "--title", help="The title to use for the plot", default=default_title
    )

    parser.add_argument(
        "--use-groups",
        help="If this flag is given, the the ORFs " "will be grouped",
        action="store_true",
    )

    args = parser.parse_args()

    orfs = bed_utils.read_bed(args.orfs)

    if args.use_groups:
        orf_lengths = [
            get_orf_lengths(orfs, ribo_utils.orf_type_labels_mapping[label])
            for label in ribo_utils.orf_type_labels
        ]

        prediction_labels = [
            latex.get_latex_safe_string(l) for l in ribo_utils.orf_type_labels
        ]

        prediction_lengths_list = orf_lengths
    else:
        orf_lengths = [
            get_orf_lengths(orfs, [orf_type]) for orf_type in ribo_utils.orf_types
        ]

        prediction_labels = [
            latex.get_latex_safe_string(l) for l in ribo_utils.orf_types
        ]

        prediction_lengths_list = orf_lengths

    if os.path.exists(args.uniprot):
        truth_nt_lengths = bio.get_uniprot_nt_lengths(args.uniprot)
        truth_label = args.uniprot_label
    else:
        truth_nt_lengths = None
        truth_label = None

    # prediction_lengths_list = [bf_lengths, chisq_lengths]
    # prediction_labels = ['BF', r'$\chi^2$']

    # input: truth_nt_lengths (array-like)
    #        prediction_lengths_list (list of array-likes)
    #        truth_label (string)
    #        prediction_labels (list of array-likes)
    #
    # if truth_nt_lengths is not defined, then the KL-divergence calculations
    # will be skipped (and it will not be shown)

    fontsize = 20
    legend_fontsize = 20
    title_fontsize = 20
    linewidth = 4

    # plot the empirical distribution of ORF lengths
    hist_min = 200
    hist_max = 5250
    hist_step = 200
    hist_range = (hist_min, hist_max)
    hist_bins = np.arange(hist_min, hist_max, hist_step)

    if truth_nt_lengths is not None:
        truth_hist, _ = np.histogram(
            truth_nt_lengths, bins=hist_bins, range=hist_range, density=True
        )
    else:
        truth_hist = None

    prediction_hists = []
    for prediction_lengths in prediction_lengths_list:
        prediction_hist, _ = np.histogram(
            prediction_lengths, bins=hist_bins, range=hist_range, density=True
        )
        prediction_hists.append(prediction_hist)

    # now, normalize the histograms
    if truth_hist is not None:
        truth_hist = truth_hist / np.sum(truth_hist)
        truth_hist += 1e-3

    for i, prediction_hist in enumerate(prediction_hists):
        prediction_hists[i] = prediction_hist / np.sum(prediction_hist)
        prediction_hists[i] += 1e-3

    kls = []
    if truth_hist is not None:
        for i, prediction_hist in enumerate(prediction_hists):
            kl = math_utils.calculate_symmetric_kl_divergence(
                truth_hist, prediction_hist, scipy.stats.entropy
            )
            kls.append(kl)

            # and update the label
            prediction_labels[i] = "{}, KL: ${:.2f}$".format(prediction_labels[i], kl)

    if truth_hist is not None:
        truth_hist = 100 * truth_hist

    for i, prediction_hist in enumerate(prediction_hists):
        prediction_hists[i] *= 100

    fig, ax = plt.subplots(figsize=(10, 5))

    cm = plt.cm.gist_earth

    x = np.arange(len(hist_bins) - 1)

    truth_cm_offset = 0.1
    if truth_hist is not None:
        color = cm(truth_cm_offset)
        ax.plot(x, truth_hist, label=truth_label, linewidth=linewidth, color=color)

    color_range = 1 - 2 * truth_cm_offset
    for i, prediction_hist in enumerate(prediction_hists):
        color = i / len(prediction_hists) * color_range
        color += 2 * truth_cm_offset
        color = cm(color)
        ax.plot(
            x,
            prediction_hist,
            label=prediction_labels[i],
            linewidth=linewidth,
            color=color,
        )

    ax.set_xlabel("Length (bp)", fontsize=fontsize)
    ax.set_ylabel("\% of predicted ORFs", fontsize=fontsize)

    if len(args.title) > 0:
        ax.set_title(args.title, fontsize=fontsize)

    ax.set_xticks(x[::2])
    ax.set_xticklabels(hist_bins[::2], fontsize=fontsize, rotation=90)

    ax.set_ylim((0, 20))
    ax.set_xlim((0, len(hist_bins)))

    # hide the "0" tick label
    yticks = ax.yaxis.get_major_ticks()
    yticks[0].label1.set_visible(False)

    # chop off everything from 3000 on
    index_of_3000 = 14
    ax.set_xlim((0, index_of_3000))
    # ax.set_xlim((0, len(uniprot_hist)-1))

    lgd = ax.legend(
        loc="center right", fontsize=legend_fontsize, bbox_to_anchor=(1.75, 0.5)
    )
    ax.tick_params(axis="both", which="major", labelsize=fontsize)

    fig.savefig(args.out, bbox_inches="tight", bbox_extra_artists=(lgd,))


if __name__ == "__main__":
    main()
