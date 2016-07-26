#! /usr/bin/env python3

import argparse
import logging
import numpy as np
import pandas as pd

import misc.bio as bio
import misc.parallel as parallel
import misc.utils as utils

import riboutils.ribo_utils as ribo_utils

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

default_seqname_prefix = ""


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script performs a permutation test to estimate the "
        "significance of the number of exact matches of Mackowiak ORFs found "
        "among the predicted ORFs and those among a random selection of ORFs."
        "Specifically, though, this script writes the number of exact matches "
        "found during each permutation to a file.\n\nThis script has two modes "
        "of operation for selecting the \"background\" for the permutation "
        "test. In the first, all ORFs are used in the background; in the "
        "second, only those ORFs which pass the basic filtering criteria (min "
        "length, min profile, more reads in the first reading frame for this "
        "sample) are used; in the second, all ORFs are considered in the "
        "background.")

    parser.add_argument('bf', help="The Bayes factor file from "
        "estimate-orfs-bayes-factors")
    parser.add_argument('mackowiak_orfs', help="The ORFs given BED12 format. This should "
        "be the output from the mackowiak2bed.py script (or that output subsequently "
        "processed); it should *not* be the raw supplementary files.")
    parser.add_argument('out', help="The output file, which is a csv file that "
        "contains the number of exact matches from each permutation")

    parser.add_argument('-a', '--use-all-orfs', help="If this flag is given, then all "
        "ORFs will be used as the background. Otherwise, only those which pass the "
        "basic filtering criteria will be used.", action='store_true')

    parser.add_argument('--seqname-prefix', help="If present, this string is prepended "
        "to all of the ORF seqnames before matching to the Mackowiak ORFs.", 
        default=default_seqname_prefix)


    parser.add_argument('-n', '--num-random-samples', help="The number of random samples "
        "to draw for the permuation test.", type=int, default=default_num_random_samples)
    
    parser.add_argument('-p', '--num-cpus', help="The number of processors to use "
        "for parallelization the sampling", type=int, default=default_num_cpus)

    parser.add_argument('-g', '--num-groups', help="The number of groups to use for "
        "performing the permutation tests. More groups means the progress bar is "
        "updated more frequently but incurs more overhead because of the parallel "
        "calls.", type=int, default=default_num_groups)    

    parser.add_argument('--seed', help="The random seed", type=int, default=default_seed)

    parser.add_argument('--min-length', help="The minimum length to predict an ORF "
        "as translated", type=int, default=default_min_length)

    parser.add_argument('--min-profile', help="ORFs with profile sum (i.e., number "
        "of reads) less than this value will not be processed.", type=float, 
        default=default_min_profile)

    
    parser.add_argument('--min-bf-mean', help="The minimum Bayes' factor mean to predict "
        "an ORF as translated (use --help for more details)", 
        type=float, default=default_min_bf_mean)
    parser.add_argument('--max-bf-var', help="The maximum Bayes' factor variance to predict "
        "an ORF as translated (use --help for more details)", 
        type=float, default=default_max_bf_var)

    parser.add_argument('--min-bf-likelihood', help="If given, then this is taken a threshold "
        "on the likelihood of translation (use --help for more details)", 
        type=float, default=default_min_bf_likelihood)

    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    np.random.seed(args.seed)

    msg = "Reading the Bayes factor file"
    logger.info(msg)

    bf = bio.read_bed(args.bf)

    msg = "Getting the filters"
    logger.info(msg)

    m_base_filter = ribo_utils.get_base_filter(bf, args.min_profile, 
                                        args.min_length)

    m_bf_filter = ribo_utils.get_bf_filter(bf, 
                                        min_bf_mean=args.min_bf_mean,
                                        max_bf_var=args.max_bf_var, 
                                        min_bf_likelihood=args.min_bf_likelihood)

    # we will choose the same number of ORFs as predicted
    m_predicted = m_base_filter & m_bf_filter
    num_predicted_orfs = sum(m_predicted)

    if args.use_all_orfs:
        bf_background = bf
    else:
        bf_background = bf[m_base_filter]

    msg = "Preparing for parallel calls"
    logger.info(msg)

    # how many samples per group?
    samples_per_group = int(np.ceil(args.num_random_samples / args.num_groups))

    # we need a dummy iterator for each group
    dummy_iter = np.arange(args.num_groups)

    samples_l = parallel.apply_parallel_iter(
        dummy_iter,
        args.num_cpus,
        ribo_utils.get_mackowiak_background,
        samples_per_group, bf_background, num_predicted_orfs,
        args.mackowiak_orfs, args.seqname_prefix,
        progress_bar=True, num_groups=args.num_groups)

    all_samples = np.concatenate(samples_l)

    samples_df = pd.DataFrame()
    samples_df['mackowiak_exact_matches'] = all_samples
    samples_df['sample_index'] = np.arange(len(all_samples))

    utils.write_df(samples_df, args.out, index=False)

if __name__ == '__main__':
    main()
