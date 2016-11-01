#! /usr/bin/env python3

import argparse
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import misc.mpl_utils as mpl_utils
import misc.np_utils as np_utils

import pickle

import logging
import misc.logging_utils as logging_utils
logger = logging.getLogger(__name__)

default_title = ""
default_min_weight = 0.001

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script visualizes the clusters found with "
        "cluster-subcodon-counts.")

    parser.add_argument('pkl', help="The pickled model file created by "
        "cluster-subcodon-counts")
    parser.add_argument('out', help="The output image")

    parser.add_argument('--title', help="The title for the plot", 
        default=default_title)
    parser.add_argument('--min-weight', help="The minimum weight required to "
        "show the associated cluster", type=float, default=default_min_weight)
    parser.add_argument('--log', help="If this flag is given, then the plot "
        "will use a log scale", action='store_true')
    
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "Reading model pickle file"
    logger.info(msg)
    model_pkl = pickle.load(open(args.pkl, 'rb'))

    msg = "Extracting clusters with minimum weight"
    logger.info(msg)

    it = enumerate(zip(model_pkl[0], model_pkl[1]))

    periodic_clusters = []

    total_weight = 0
    for i, (m, w) in it:
        if w > args.min_weight:        
            total_weight += w
            periodic_clusters.append(i)

    msg = "Finding linear best fit line"
    logger.info(msg)

    c = model_pkl[0][periodic_clusters, 0]
    x = model_pkl[0][periodic_clusters, 1]
    y = model_pkl[0][periodic_clusters, 2]
    s = model_pkl[1][periodic_clusters]

    fit = np_utils.fit_with_least_squares(x, y, w=s)
    (slope, intercept, power, r_sqr) = fit

    msg = "Plotting clusters"
    logger.info(msg)

    min_val = min(min(x), min(y)) * 0.8
    max_val = max(max(x), max(y)) * 1.2
    lim = (min_val, max_val)

    fig, ax = plt.subplots()

    # axes and labels and things
    ax.set_aspect('equal')
    ax.set_xlabel("Frame +1")
    ax.set_ylabel("Frame +2")

    ax.set_xlim(lim)
    ax.set_ylim(lim)

    if args.log:
        ax.set_xscale('log')
        ax.set_yscale('log')

    cm = plt.cm.Blues

    norm = None
    if args.log:
        norm=matplotlib.colors.LogNorm()
        
    sc = ax.scatter(x, y, c=c, cmap=cm, s=s*1000, norm=norm)
    cb = plt.colorbar(sc, ax=ax)
    cb.set_label("In-frame")

    text = "Accounts for {:.0%} of likelihood".format(total_weight)
    ax.annotate(text, (0.25,0.75), xycoords='axes fraction')

    # draw the fit line
    mpl_utils.plot_trend_line(ax, x, intercept, slope, power)

    # write the fit information
    rsqr_str = "$R^2$ = {:.2f}".format(r_sqr)
    slope_str = "slope = {:.2f}".format(slope)
    intercept_str = "intercept = {:.2f}".format(intercept)
    strs = [rsqr_str, slope_str, intercept_str]
    text = '\n'.join(strs)

    ax.annotate(text, (0.55, 0.15), xycoords='axes fraction')

    if len(args.title) > 0:
        ax.set_title(args.title)
        
    msg = "Writing the plot to disk"
    logger.info(msg)
        
    fig.savefig(args.out)

if __name__ == '__main__':
    main()
