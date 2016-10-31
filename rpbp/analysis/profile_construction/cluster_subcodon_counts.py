#! /usr/bin/env python3

import matplotlib
matplotlib.use('agg')

import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
import scipy
import sklearn.mixture

import misc.bio_utils.bed_utils as bed_utils
import misc.logging_utils as logging_utils
import misc.np_utils as np_utils

import logging
import misc.logging_utils as logging_utils
logger = logging.getLogger(__name__)

default_max_iter=100
default_n_components = 100
default_seed=8675309
default_min_weight = 0.01

default_fraction = 0.1

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script clusters the ORFs based on their subcodon "
        "counts using a DP-GMM; the means and weights of the clusters are "
        "written to a pickle file which consists of a list. The first element "
        "of the list are the means, and the second is the weights.")

    parser.add_argument('bf', help="The bayes factor file containing counts")
    parser.add_argument('out', help="The output (pickle) file")

    parser.add_argument('--fraction', help="The top <fraction> genes, based "
        "on normalized read counts, will be used for clustering", type=float,
        default=default_fraction)

    parser.add_argument('--max-iter', help="The maximum number of iterations "
        "for clustering", type=int, default=default_max_iter)
    parser.add_argument('--n-components', help="The maximum number of "
        "clusters", type=int, default=default_n_components)
    parser.add_argument('--seed', help="The seed for the random number "
        "generator", type=int, default=default_seed)
    
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "Reading BF file"
    logger.info(msg)
    bf = bed_utils.read_bed(args.bf)

    msg = "Extracting top k% of ORFs"
    logger.info(msg)
    
    # calculate the normalized read coverage
    total_read_coverage = bf["x_1_sum"] + bf["x_2_sum"] + bf["x_3_sum"]
    rpk = total_read_coverage / bf['orf_len']
    sorted_rpk_indices = np.argsort(rpk)

    # and get the best args.fraction of them
    num_orfs = int(len(rpk) * args.fraction)
    top_k_orfs = sorted_rpk_indices.tail(num_orfs)

    msg = "Finding subcodon clusters"
    logger.info(msg)

    x_i_fields = ["x_1_sum", "x_2_sum", "x_3_sum"]
    X = bf.iloc[top_k_orfs.values][x_i_fields]

    model = np_utils.fit_bayesian_gaussian_mixture(X, 
        max_iter=args.max_iter, n_components=args.n_components, seed=args.seed)

    msg = "Writing means and weights to disk"
    logger.info(msg)

    to_pkl = [model.means_, model.weights_]
    pickle.dump(to_pkl, open(args.out, 'wb'))

if __name__ == '__main__':
    main()
