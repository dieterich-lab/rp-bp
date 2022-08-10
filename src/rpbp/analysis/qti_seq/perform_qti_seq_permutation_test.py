#! /usr/bin/env python3

import argparse
import logging

import pbiotools.misc.utils as utils

logger = logging.getLogger(__name__)

default_num_cpus = 1
default_num_groups = 500
default_num_random_samples = 10000

default_seed = 8675309

default_min_bf_mean = 5
default_min_bf_likelihood = 0.5
default_max_bf_var = None
default_min_length = 20
default_min_profile = 5

default_max_length = 300

default_tmp = None


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script performs a permutation test to estimate the"
        "significance of the number of micropeptides which have a nearby QTI-"
        "seq peak. Specifically, this script writes the number of micropeptides "
        "with a nearby QTI-seq peak found during each permutation to a file.",
    )

    parser.add_argument(
        "bf", help="The Bayes factor file from " "estimate-orfs-bayes-factors"
    )
    parser.add_argument(
        "qti_seq_peaks",
        help="The QTI-seq peaks in BED format. "
        "This should be the output of qtipeaks2bed.py, not the Supplemental "
        "Files distributed with the paper",
    )
    parser.add_argument(
        "out",
        help="The output (csv.gz) file. Each line will "
        "contain a permutation number and the number of micropeptides with a "
        "QTI-seq peak.",
    )

    parser.add_argument(
        "--max-length",
        help="The maximum length to consider an "
        'ORF as translated. This is used to keep only "micropeptides".',
        type=int,
        default=default_max_length,
    )

    parser.add_argument(
        "-n",
        "--num-random-samples",
        help="The number of random samples " "to draw for the permuation test.",
        type=int,
        default=default_num_random_samples,
    )

    parser.add_argument(
        "-p",
        "--num-cpus",
        help="The number of processors to use " "for parallelization the sampling",
        type=int,
        default=default_num_cpus,
    )

    parser.add_argument(
        "-g",
        "--num-groups",
        help="The number of groups to use for "
        "performing the permutation tests. More groups means the progress bar is "
        "updated more frequently but incurs more overhead because of the parallel "
        "calls.",
        type=int,
        default=default_num_groups,
    )

    parser.add_argument(
        "--seed", help="The random seed", type=int, default=default_seed
    )

    parser.add_argument(
        "--min-length",
        help="The minimum length to predict an ORF " "as translated",
        type=int,
        default=default_min_length,
    )

    parser.add_argument(
        "--min-profile",
        help="ORFs with profile sum (i.e., number "
        "of reads) less than this value will not be processed.",
        type=float,
        default=default_min_profile,
    )

    parser.add_argument(
        "--min-bf-mean",
        help="The minimum Bayes' factor mean to predict "
        "an ORF as translated (use --help for more details)",
        type=float,
        default=default_min_bf_mean,
    )
    parser.add_argument(
        "--max-bf-var",
        help="The maximum Bayes' factor variance to predict "
        "an ORF as translated (use --help for more details)",
        type=float,
        default=default_max_bf_var,
    )

    parser.add_argument(
        "--min-bf-likelihood",
        help="If given, then this is taken a threshold "
        "on the likelihood of translation (use --help for more details)",
        type=float,
        default=default_min_bf_likelihood,
    )

    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)


if __name__ == "__main__":
    main()
