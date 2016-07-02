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
    orf_type_groups = orfs.groupby('orf_type')
    counts = orf_type_groups.size()

    if args.use_groups:
        labels = ribo_utils.orf_type_labels
        fracs = [get_orf_label_counts(counts, l) for l in labels]
    else:
        fracs = counts.values
        labels = np.array(counts.index)

    labels = ["{} ({})".format(l,f) for l,f in zip(labels, fracs)]

    fig, ax = plt.subplots(figsize=(5,5))

    cmap = plt.cm.Blues
    colors = cmap(np.linspace(0., 1., len(labels)))

    #patches, texts, autotexts = ax.pie(fracs, labels=labels, colors=colors, autopct="%.1f%%")
    patches, texts = ax.pie(fracs, colors=colors)
    lgd = ax.legend(patches, labels, loc="center right", bbox_to_anchor=(0,0.5))

    #for autotext, count in zip(autotexts, fracs):
        #autotext.set_text(u"%s (%d)" % (autotext.get_text(), counts))
        #autotext.set_text("{:.1e}".format(count))

    if len(args.title) > 0:
        ax.set_title(args.title)

    fig.savefig(args.out, bbox_extra_artists=(lgd,), bbox_inches='tight')

if __name__ == '__main__':
    main()
