#! /usr/bin/env python3

import matplotlib
matplotlib.use('agg')

import argparse

import matplotlib.pyplot as plt
import numpy as np

import misc.bio as bio

import riboutils.ribo_utils as ribo_utils


default_title = ""

def get_orf_label_counts(counts, orf_label):
    orf_types = ribo_utils.orf_type_labels_mapping[orf_label]
    orf_label_counts = np.sum(counts[orf_types])
    
    if np.isnan(orf_label_counts):
        orf_label_counts = 0
    return orf_label_counts

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script creates a pie chart which shows the proportion of "
        "each ORF type in a given BED12+ file. Optionally, the ORFs can be grouped "
        "into similar types.")

    parser.add_argument('orfs', help="The BED12+ file with the ORFs")
    parser.add_argument('out', help="The output (image) file")

    parser.add_argument('--title', help="The title to use for the plot", 
        default=default_title)

    parser.add_argument('--use-groups', help="If this flag is given, the the ORFs "
        "will be grouped", action='store_true')
    
    args = parser.parse_args()

    orfs = bio.read_bed(args.orfs)

    strands = ['+', '-']
    fracs = []
    labels = []
    for strand in ['+', '-']:
        m_strand = orfs['strand'] == strand
        orf_type_groups = orfs[m_strand].groupby('orf_type')
        counts = orf_type_groups.size()

        if args.use_groups:
            lab = ribo_utils.orf_type_labels
            fr = [get_orf_label_counts(counts, l) for l in lab]
        else:
            fr = counts.values
            lab = np.array(counts.index)

        lab = ["{} ({})".format(l,f) for l,f in zip(lab, fr)]

        fracs.append(fr)
        labels.append(lab)

    fig, axes = plt.subplots(ncols=2, figsize=(10,5))

    cmap = plt.cm.Blues
    colors = cmap(np.linspace(0., 1., len(labels[0])))

    # forward strand ORFs
    patches, texts = axes[0].pie(fracs[0], colors=colors)
    lgd_1 = axes[0].legend(patches, labels[0], loc="center right", bbox_to_anchor=(0,0.5))
    axes[0].set_title("Strand: {}".format(strands[0]))

    # reverse strand ORFs
    patches, texts = axes[1].pie(fracs[1], colors=colors)
    lgd_2 = axes[1].legend(patches, labels[1], loc="center right", bbox_to_anchor=(2.0,0.5))
    axes[1].set_title("Strand: {}".format(strands[1]))

    if len(args.title) > 0:
        fig.suptitle(args.title)

    fig.savefig(args.out, bbox_extra_artists=(lgd_1, lgd_2), bbox_inches='tight')

if __name__ == '__main__':
    main()
