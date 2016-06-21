#! /usr/bin/env python3

import argparse
import json
import pickle
import logging

import time

import numpy as np
import pandas as pd
import scipy.io

import misc.bio as bio
import misc.parallel as parallel
import misc.utils as utils

default_orf_num_field = 'orf_num'
default_orf_type_field = 'orf_type'

default_orf_types = []

default_num_orfs = 0
default_num_cpus = 1
default_num_groups = 100

default_iterations = 200
default_chains = 2
default_seed = 8675309

def get_bayes_factor(profile, translated_models, untranslated_models, args):
    """ This function calculates the Bayes' factor for a single ORF profile. 

        Args:
            profile (np.array): the (dense) profile for this ORF

            translated_models (list of pystan.StanModel): the models which explain translation

            untranslated_models (list of pystan.StanModel): the models which account for background

            args (namespace): a namespace (presumably from argparse) which includes the following:
                seed (int): random seed for initializing MCMC
                chains (int): the number of MCMC chains
                iterations (int): the number of iterations for each chain

        Returns:
            pd.Series: a series containing:
            
                the mean and variance for each of the following estimated values:
                    bayes_factor
                    p_translated
                    p_background
                    translated_location
                    translated_scale
                    background_location
                    background_scale

                the chi-square p-value
                the delta_l and delta_h values
    """
    profile_sum = sum(profile)
    
    # split the signal based on frame
    x_1 = profile[0::3]
    x_2 = profile[1::3]
    x_3 = profile[2::3]
    T = len(x_1)
    nonzero_x_1 = np.count_nonzero(x_1)
   
    x_1_sum = sum(x_1)
    x_2_sum = sum(x_2)
    x_3_sum = sum(x_3)

    ret = {
        "p_translated_mean": float('-inf'),
        "p_translated_var": float('-inf'),
        "p_background_mean": float('-inf'),
        "p_background_var": float('-inf'),
        "translated_location_mean": float('-inf'),
        "translated_location_var": float('-inf'),
        "translated_scale_mean": float('-inf'),
        "translated_scale_var": float('-inf'),
        "background_location_mean": float('-inf'),
        "background_location_var": float('-inf'),
        "background_scale_mean": float('-inf'),
        "background_scale_var": float('-inf'),
        "bayes_factor_mean": float('-inf'),
        "bayes_factor_var": float('-inf'),
        "chi_square_p": float('-inf'),
        "delta_l": float('-inf'),
        "delta_h": float('-inf'),
        "x_1_sum": x_1_sum,
        "x_2_sum": x_2_sum,
        "x_3_sum": x_3_sum,
        "profile_sum": profile_sum
    }
    ret = pd.Series(ret)

    # first, make sure this profile was not filtered during smoothing
    if profile_sum == 0:
        return ret

     
    # chi-square values
    f_obs = [x_1_sum, x_2_sum, x_3_sum]
    chisq, chi_square_p = scipy.stats.chisquare(f_obs)
    ret['chi_square_p'] = chi_square_p
 
    # check if we only wanted the chi square value
    if args.chi_square_only:
        return ret
      
    # check if something odd happens with the length
    # this should already be checked before calling the function.
    if (T != len(x_2)) or (T != len(x_3)):
        return ret

    # construct the input for Stan
    data = {
        "x_1": x_1,
        "x_2": x_2,
        "x_3": x_3,
        "T": T,
        "nonzero_x_1": nonzero_x_1
    }

    m_translated = [tm.sampling(data=data, iter=args.iterations, chains=args.chains, n_jobs=1, 
        seed=args.seed, refresh=0) for tm in translated_models]
    
    m_background = [bm.sampling(data=data, iter=args.iterations, chains=args.chains, n_jobs=1, 
        seed=args.seed, refresh=0) for bm in untranslated_models]

    
    # extract the parameters of interest
    m_translated_ex = [m.extract(pars=['lp__', 'background_location', 'background_scale', 'delta_l', 'delta_h']) 
        for m in m_translated]

    m_background_ex = [m.extract(pars=['lp__', 'background_location', 'background_scale']) 
        for m in m_background]

    # now, choose the best model of each class,  based on mean likelihood
    m_translated_means = [np.mean(m_ex['lp__']) for m_ex in m_translated_ex]
    m_background_means = [np.mean(m_ex['lp__']) for m_ex in m_background_ex]

    max_translated_mean = np.argmax(m_translated_means)
    max_background_mean = np.argmax(m_background_means)

    # select the best sampling results
    m_translated_ex = m_translated_ex[max_translated_mean]
    m_background_ex = m_background_ex[max_background_mean]

    # extract the relevant means and variances
    ret['p_translated_mean'] = np.mean(m_translated_ex['lp__'])
    ret['p_translated_var'] = np.var(m_translated_ex['lp__'])

    ret['p_background_mean'] = np.mean(m_background_ex['lp__'])
    ret['p_background_var'] = np.var(m_background_ex['lp__'])

    ret['translated_location_mean'] = np.mean(m_translated_ex['background_location'])
    ret['translated_location_var'] = np.var(m_translated_ex['background_location'])

    ret['translated_scale_mean'] = np.mean(m_translated_ex['background_scale'])
    ret['translated_scale_var'] = np.var(m_translated_ex['background_scale'])

    ret['background_location_mean'] = np.mean(m_background_ex['background_location'])
    ret['background_location_var'] = np.var(m_background_ex['background_location'])

    ret['background_scale_mean'] = np.mean(m_background_ex['background_scale'])
    ret['background_scale_var'] = np.var(m_background_ex['background_scale'])

    ret['delta_l'] = np.mean(m_translated_ex['delta_l'])
    ret['delta_h'] = np.mean(m_translated_ex['delta_h'])

    # the (log of) the Bayes factor is the difference between two normals:
    # (the best translated model) - (the best background model)
    #
    # thus, it is also a normal whose mean is the difference of the two means
    # and whose variance is the sum of the two variances
    ret['bayes_factor_mean'] = ret['p_translated_mean'] - ret['p_background_mean']
    ret['bayes_factor_var'] = ret['p_translated_var'] + ret['p_background_var']
    
    return ret

def get_all_bayes_factors(orfs, args):
    """ This function calculates the Bayes' factor term for each region in regions. See the
        description of the script for the Bayes' factor calculations.

        Args:
            orfs (pd.DataFrame) : a set of orfs. The columns must include:
                orf_num
                exon_lengths

            args (namespace) : a namespace containing the models and profiles filenames
            
        Returns:
            pandas.Series: the Bayes' factors (and other estimated quantities) for each region
    """

    # read in the signals and sequences
    logging.debug("Reading profiles")
    profiles = scipy.io.mmread(args.profiles).tocsr()
    
    logging.debug("Reading models")
    translated_models = [pickle.load(open(tm, 'rb')) for tm in args.translated_models]
    untranslated_models = [pickle.load(open(bm, 'rb')) for bm in args.untranslated_models]

    logging.debug("Applying on regions")
    bfs = []
    for idx, row in orfs.iterrows():
        orf_num = row[args.orf_num_field]
        orf_len = row['orf_len']

        # sometimes the orf_len is off...
        if orf_len % 3 != 0:
            msg = "Found an ORF whose length was not 0 mod 3. Skipping. orf_id: {}".format(row['id'])
            logging.warn(msg)
            continue

        profile = utils.to_dense(profiles, orf_num, float, length=orf_len)

        row_bf = get_bayes_factor(profile, translated_models, untranslated_models, args)
        row = row.append(row_bf)

        bfs.append(row)

    bfs = pd.DataFrame(bfs)
    return bfs

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""
            This script uses Hamiltonian MCMC with Stan to estimate translation parameters
            for a set of regions (presumably ORFs). Roughly, it takes as input:

            (1) a set of regions (ORFs) and their corresponding profiles
            (2) a "translated" model which gives the probability that a region is translated
            (3) an "untranslated" model which gives the probability that a region is not translated

            The script calculates both the Bayes' factor and \chi^2 value for each ORF.
        """
        )

    parser.add_argument('profiles', help="The ORF profiles (counts) (mtx)")
    parser.add_argument('regions', help="The regions (ORFs) for which predictions will "
        "be made (BED12+)")
    
    parser.add_argument('out', help="The output file for the Bayes' factors (BED12+)")

    parser.add_argument('--chi-square-only', help="If this flag is present, then only the chi "
        "square test will be performed for each ORF. This can also be a way to get the counts "
        "within each of the ORFs.", action='store_true')
    
    parser.add_argument('--translated-models', help="The models to use as H_t (pkl)", nargs='+')
    parser.add_argument('--untranslated-models', help="The models to use as H_u (pkl)", nargs='+')
    
    parser.add_argument('--orf-types', help="If values are given, then only orfs with "
        "those types are processed.", nargs='*', default=default_orf_types)
    parser.add_argument('--orf-type-field', default=default_orf_type_field)

    parser.add_argument('-s', '--seed', help="The random seeds to use for inference",
        type=int, default=default_seed)
    parser.add_argument('-c', '--chains', help="The number of MCMC chains to use", type=int,
        default=default_chains)
    parser.add_argument('-i', '--iterations', help="The number of MCMC iterations to use for "
        "each chain", type=int, default=default_iterations)
    
    parser.add_argument('--num-orfs', help="If n>0, then only this many ORFs will be processed",
        type=int, default=default_num_orfs)
    parser.add_argument('--orf-num-field', default=default_orf_num_field)

    parser.add_argument('--num-cpus', help="The number of CPUs to use. ", type=int,
        default=default_num_cpus)
    parser.add_argument('--do-not-compress', help="Unless otherwise specified, the output will "
        "be written in GZip format", action='store_true')

    parser.add_argument('-g', '--num-groups', help="The number of groups into which to split "
        "the ORFs. More groups means the progress bar is updated more frequently but incurs "
        "more overhead because of the parallel calls.", type=int, default=default_num_groups)

    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    # read in the regions and apply the filters
    msg = "Reading and filtering ORFs"
    logging.info(msg)
    regions = bio.read_bed(args.regions)
    
    if args.num_orfs > 0:
        mask_orfs = regions[args.orf_num_field] < args.num_orfs
        regions = regions[mask_orfs]

    # add the lengths of the orfs
    exon_lengths = regions.apply(bio.get_bed_12_feature_length, axis=1)
    regions['orf_len'] = exon_lengths

    if len(args.orf_types) > 0:
        mask_a = regions[args.orf_type_field].isin(args.orf_types)
        regions = regions[mask_a]

    regions = regions.reset_index(drop=True)

    msg = "Number of regions after filtering: {}".format(len(regions))
    logging.info(msg)

    # read in everything else in the parallel call

    # calculate the bayes' factor for each region
    bfs_l = parallel.apply_parallel_split(regions, args.num_cpus, get_all_bayes_factors, args,
        num_groups=args.num_groups, progress_bar=True)
    bfs = pd.concat(bfs_l)

    # write the results as a bed12+ file
    bio.write_bed(bfs, args.out)

if __name__ == '__main__':
    main()

