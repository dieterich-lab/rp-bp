#! /usr/bin/env python3

import matplotlib
matplotlib.use('agg')

import argparse
import logging
import os
import yaml

import matplotlib.pyplot as plt
import numpy as np
import misc.bio as bio
import misc.utils as utils

import riboutils.ribo_filenames as filenames
import riboutils.ribo_utils as ribo_utils

logger = logging.getLogger(__name__)

default_min_rpkm = 0
default_max_rpkm = 5

default_min_bf = -10000
default_max_bf = 1000

default_title = ""

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script plots the (log) Bayes factor against the estimated "
        "RPKM for all ORFs. All relevant values will be clipped according to the "
        "specified arguments for viewing.")

    parser.add_argument('config', help="The (yaml) config file")
    parser.add_argument('name', help="The name of the dataset or replicate to plot")
    parser.add_argument('out', help="The output image file")

    parser.add_argument('-p', '--use-predictions', help="If this flag is present, then "
        "the \"predicted ORFs\" files will be used. Otherwise, all ORFs in the dataset "
        "will be visualized.", action='store_true')
    parser.add_argument('-r', '--is-replicate', help="If the name corresponds to one "
        "of the replicates, this flag must be used to ensure the filenames are "
        "handled correctly.", action='store_true')

    parser.add_argument('--title', default=default_title)

    parser.add_argument('--min-rpkm', type=float, default=default_min_rpkm)
    parser.add_argument('--max-rpkm', type=float, default=default_max_rpkm)
    parser.add_argument('--min-bf', type=float, default=default_min_bf)
    parser.add_argument('--max-bf', type=float, default=default_max_bf)
    
    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    config = yaml.load(open(args.config))
    note = config.get('note', None)

    if args.is_replicate:
        lengths = None
        offsets = None
    else:
        lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(config, args.name)
        
    fraction = config.get('smoothing_fraction', None)
    reweighting_iterations = config.get('smoothing_reweighting_iterations', None)

    # we will need these to get the appropriate log BFs
    if args.use_predictions:
        bayes_factors = filenames.get_riboseq_predicted_orfs(config['riboseq_data'], args.name,
            length=lengths, offset=offsets, is_unique=True, note=note, is_smooth=True,
            fraction=fraction, reweighting_iterations=reweighting_iterations)
    else:
        bayes_factors = filenames.get_riboseq_bayes_factors(config['riboseq_data'], args.name,
            length=lengths, offset=offsets, is_unique=True, note=note, is_smooth=True,
            fraction=fraction, reweighting_iterations=reweighting_iterations)

    if not os.path.exists(bayes_factors):
        msg = ("Could not find the Bayes factor file: {}\nIf this is for a particular "
            "sample and the --merge-replicates option was used, this is not a problem. "
            "Will not create this scatter plot".format(bayes_factors))
        logger.warning(msg)
        return

    msg = "Reading Bayes factors"
    logger.info(msg)
    bayes_factors = bio.read_bed(bayes_factors)

    # we need these to get the raw counts for calculating RPKM

    # we always need all of the counts, so no need to check which ORFs
    rpchi_pvalues = filenames.get_riboseq_bayes_factors(config['riboseq_data'], args.name, 
        length=lengths, offset=offsets, is_unique=True, note=note, is_smooth=False)

    
    if not os.path.exists(rpchi_pvalues):
        msg = ("Could not find the Rp-chi pvalues file: {}\nIf this is for a particular "
            "sample and the --merge-replicates option was used, this is not a problem. "
            "Will not create this scatter plot".format(rpchi_pvalues))
        logger.warning(msg)
        return

    msg = "Reading Rp-chi pvalues"
    logger.info(msg)
    rpchi_pvalues = bio.read_bed(rpchi_pvalues)

    msg = "Calculating RPKM values"
    logger.info(msg)

    # we approximate the number of mapping reads as the sum across all ORFs.
    # this double-counts some reads
    num_reads = np.sum(rpchi_pvalues['profile_sum'])
    all_rpkm = (1e6 * rpchi_pvalues['x_1_sum']) / (rpchi_pvalues['orf_len'] * num_reads)

    # only include things that have some reads in the visualization
    m_rpkm = all_rpkm > 0

    msg = "Creating plot"
    logger.info(msg)

    fig, ax = plt.subplots(figsize=(10, 5))

    cm = plt.cm.gist_earth

    for i, orf_label in enumerate(ribo_utils.orf_type_labels):
        
        orf_types = ribo_utils.orf_type_labels_mapping[orf_label]
        m_type = bayes_factors['orf_type'].isin(orf_types)

        # now, pull out the RPKMs
        if args.use_predictions:
            # if we are using predictions, we have to filter and join
            orf_ids = bayes_factors.loc[m_rpkm & m_type, 'id']
            bfs = np.array(bayes_factors.loc[m_rpkm & m_type, 'bayes_factor_mean'])

            m_ids = rpchi_pvalues['id'].isin(orf_ids)
            rpkm = np.array(all_rpkm[m_ids])

        else:
            # otherwise ,the data frames match, so we can just use the masks
            rpkm = np.array(all_rpkm[m_rpkm & m_type])
            bfs = np.array(bayes_factors.loc[m_rpkm & m_type, 'bayes_factor_mean'])
        
        rpkm = np.clip(rpkm, args.min_rpkm, args.max_rpkm)
        bfs = np.clip(bfs, args.min_bf, args.max_bf)
        
        color = i / len(ribo_utils.orf_type_labels)
        color = cm(color)
        
        label = "{} ({})".format(orf_label, len(rpkm))

        ax.scatter(rpkm, bfs, label=label, color=color, edgecolor='k')


    ax.set_ylim((args.min_bf * 1.5, args.max_bf * 1.5))
    ax.set_xlim((args.min_rpkm * 1.5, args.max_rpkm * 1.25))

    ax.set_yscale('symlog')
    ax.set_xscale('symlog')

    ax.set_xlabel('RPKM')
    ax.set_ylabel('log BF')

    lgd = ax.legend(loc='center right', bbox_to_anchor=(1.5, 0.5))

    if len(args.title) > 0:
        ax.set_title(args.title)

    fig.savefig(args.out, bbox_inches='tight', bbox_extra_artists=(lgd,))



if __name__ == '__main__':
    main()
