#! /usr/bin/env python3

import argparse
import logging
import numpy as np
import scipy.io
import scipy.sparse


import statsmodels.api as sm
lowess = sm.nonparametric.lowess

import misc.bio as bio
import misc.external_sparse_matrix_list as external_sparse_matrix_list
import misc.parallel as parallel
import misc.utils as utils

default_min_length = 0
default_max_length = 0
default_min_signal = 5

default_fraction = 0.5
default_reweighting_iterations = 0

default_num_cpus = 1
default_num_orfs = 0

def smooth_profile(profile, args):

    # first, convert it to a dense matrix
    # this is always a single row, so always convert row 0

    # HACK ALERT: finding the length of the ORF
    # this is similar in spirit to this stack overflow approach:
    # http://stackoverflow.com/questions/24841271/finding-maximum-value-and-their-indices-in-a-sparse-lil-matrix-scipy-python

    #by design, the min value of the matrix is -1, and that is at the end of the ORF

    # find the data position of -1
    pos = profile.data.argmin()
    
    # now, find the actual index of that data position
    neg_index = profile.indices[pos]

    profile = utils.to_dense(profile, 0, float, length=neg_index)
    smoothed_profile = np.zeros_like(profile)

    msg = "Length of profile: {}".format(len(profile))
    logging.debug(msg)

    # make sure the length is okay
    if (args.min_length > 0) and (neg_index < args.min_length):
        return scipy.sparse.csr_matrix(smoothed_profile)
    
    if (args.max_length > 0) and (neg_index > args.max_length):
        return scipy.sparse.csr_matrix(smoothed_profile)

    # split the signal based on frame
    x_1 = profile[0::3]
    x_2 = profile[1::3]
    x_3 = profile[2::3]

    # do not bother smoothing signals for which we will not estimate the Bayes factor
    # in these cases, just return a 0 signal
    x_1_sum = sum(x_1)
    x_2_sum = sum(x_2)
    x_3_sum = sum(x_3)
    
    # make sure there is enough signal to process
    if x_1_sum < args.min_signal:
        return scipy.sparse.csr_matrix(smoothed_profile)

    # also, make sure we have more reads in x_1 than the others
    if x_1_sum < (x_2_sum + x_3_sum):
        return scipy.sparse.csr_matrix(smoothed_profile)


    exog = np.arange(len(x_1))

    # x_1
    endog = x_1
    smoothed_x_1 = lowess(endog, exog, is_sorted=True, return_sorted=False, 
        it=args.reweighting_iterations, frac=args.fraction)
    
    # x_2
    endog = x_2
    smoothed_x_2 = lowess(endog, exog, is_sorted=True, return_sorted=False, 
        it=args.reweighting_iterations, frac=args.fraction)
    
    # x_3
    endog = x_3
    smoothed_x_3 = lowess(endog, exog, is_sorted=True, return_sorted=False, 
        it=args.reweighting_iterations, frac=args.fraction)
    
    
    smoothed_profile[0::3] = smoothed_x_1
    smoothed_profile[1::3] = smoothed_x_2
    smoothed_profile[2::3] = smoothed_x_3

    sparse_smoothed_profile = scipy.sparse.csr_matrix(smoothed_profile)
    
    return sparse_smoothed_profile


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script smoothes the extracted ORF profiles. In particular, it "
        "uses LOWESS smoothing, with the specified parameters, on the counts in each "
        "frame, independently. Finally, it interleaves the smoothed frames to create the "
        "final ORF profile.\n\nPlease see statsmodels.api.nonparametric.lowess for more "
        "details about the LOWESS smoothing.")

    parser.add_argument('orfs', help="The (BED12+) ORFs file")
    parser.add_argument('profiles', help="The (mtx) ORF profiles")
    parser.add_argument('out', help="The output (mtx) file")

    parser.add_argument('--min-length', help="ORFs with length less than this value will not "
        "be processed", type=int, default=default_min_length)
    parser.add_argument('--max-length', help="ORFs with length greater than this value will not "
        "be processed", type=int, default=default_max_length)
    parser.add_argument('--min-signal', help="ORFs with profile signal less than this value "
        "will not be processed.", type=float, default=default_min_signal)


    parser.add_argument('--fraction', help="The fraction of signal to use in LOWESS", 
        type=float, default=default_fraction)

    parser.add_argument('--reweighting-iterations', help="The number of reweighting "
        "iterations to use in LOWESS. Please see the statsmodels documentation for a "
        "detailed description of this parameter.", type=int, default=default_reweighting_iterations)
     
    parser.add_argument('--num-cpus', help="The number of processes to use", 
        type=int, default=default_num_cpus)
    parser.add_argument('-k', '--num-orfs', help="If  n>0, then only the first n orfs "
        "will be processed.", type=int, default=default_num_orfs)

    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    # read in the regions and apply the filters
    msg = "Reading ORFs"
    logging.info(msg)

    orfs = bio.read_bed(args.orfs)
    
    if args.num_orfs > 0:
        mask_orfs = orfs['orf_num'] < args.num_orfs
        orfs = orfs[mask_orfs]

    orfs = orfs.sort_values('orf_num')
    
    # add the lengths of the orfs
    msg = "Getting ORF lengths"
    logging.info(msg)

    orf_lengths = parallel.apply_parallel(orfs, args.num_cpus, 
        bio.get_bed_12_feature_length, progress_bar=True)
    orf_lengths = np.array(orf_lengths)

    msg = "Reading profiles"
    logging.info(msg)

    profiles = scipy.io.mmread(args.profiles).tolil()

    if args.num_orfs > 0:
        profiles = profiles[:args.num_orfs]

        msg = "Only using first {} ORF profiles".format(args.num_orfs)
        logging.debug(msg)

    msg = "Adding end-of-ORF marker to all ORFs"
    logging.info(msg)

    for i in range(profiles.shape[0]):
        length = orf_lengths[i]
        profiles[i, length] = -1

    profiles = profiles.tocsr()
    
    msg = "Smoothing profiles"
    logging.info(msg)

    total = profiles.shape[0]
    msg = "Number of sparse row vectors: {}".format(total)
    logging.debug(msg)

    smoothed_profiles_list = parallel.apply_parallel_iter(
        profiles,  args.num_cpus, 
        smooth_profile, args,
        progress_bar=True, total=total)

    msg = "Converting smoothed profiles into sparse matrix"
    logging.info(msg)

    smoothed_profiles = external_sparse_matrix_list.to_sparse_matrix(smoothed_profiles_list)

    msg = "Writing sparse matrix to disk"
    logging.info(msg)

    scipy.io.mmwrite(args.out, smoothed_profiles)

if __name__ == '__main__':
    main()
