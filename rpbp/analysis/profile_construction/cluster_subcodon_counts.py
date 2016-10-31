#! /usr/bin/env python3

import argparse

import logging
import misc.logging_utils as logging_utils
logger = logging.getLogger(__name__)

default_max_iter=100
default_n_components = 100
default_seed=8675309
default_min_weight = 0.01

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script clusters the ORFs based on their subcodon "
        "counts. Only clusters which exhibit the basic periodicity pattern in "
        "which we are interested are kept, i.e., the mean of the in-frame "
        "counts is larger than either of the others. It visualizes the results "
        "as a scatter plot.")

    parser.add_argument('bf', help="The bayes factor file containing counts")
    parser.add_argument('out', help="The output (image) file")

    parser.add_argument('--title', help="The title for the image")
    parser.add_argument('--max-iter', help="The maximum number of iterations "
        "for clustering", type=int, default=default_max_iter)
    parser.add_argument('--n-components', help="The maximum number of "
        "clusters", type=int, default=default_n_components)
    parser.add_argument('--seed', help="The seed for the random number "
        "generator", type=int, default=default_seed)
    parser.add_argument('--min-weight', help="The minumum weight assigned to "
        "a cluster to include it in the plot", type=float, 
        default=default_min_weight)

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "Reading BF file"
    logger.info(msg)
    bf = bed_utils.read_bed(args.bf)

    msg = "Finding subcodon clusters"
    logger.info(msg)

    x_i_fields = ["x_1_sum", "x_2_sum", "x_3_sum"]
    X = bf[x_i_fields]

    model = np_utils.fit_bayesian_gaussian_mixture(X, 
        max_iter=args.max_iter, n_components=args.n_components, seed=args.seed)

    msg = "Finding the periodic clusters"
    logger.info(msg)
    it = enumerate(zip(model.means_, model.weights_))

    periodic_clusters = []

    for i, (m, w) in it:
        if w > args.min_weight:
            if (m[0] > m[1]) and (m[0] > m[2]):
                periodic_clusters.append(i)

    c = model.means_[periodic_clusters, 0]
    x = model.means_[periodic_clusters, 1]
    y = model.means_[periodic_clusters, 2]

    msg = "Creating the plot"
    logger.info(msg)

    max_val = max(max(x), max(y))
    max_val = max_val * 1.2

    fig, ax = plt.subplots()

    ax.set_aspect('equal')

    ax.set_xlabel("Frame +1")
    ax.set_ylabel("Frame +2")

    cm = plt.cm.Blues

    sc = ax.scatter(x, y, c=c, cmap=cm, s=30)
    cb = plt.colorbar(sc, ax=ax)
    cb.set_label("In-frame")

    lim = (0, max_val)

    ax.set_xlim(lim)
    ax.set_ylim(lim)

    if len(args.title) > 0:
        ax.set_title(args.title)

    msg = "Writing the plot to disk"
    logger.info(msg)

    fig.savefig(args.out)

if __name__ == '__main__':
    main()
