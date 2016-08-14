#! /usr/bin/env python3

import argparse
import logging
import os
import pickle
import sys

import numpy as np
import pandas as pd

import misc.logging_utils as logging_utils
import misc.parallel as parallel
import misc.utils as utils

logger = logging.getLogger(__name__)

default_num_cpus = 1

default_periodic_models = []
default_nonperiodic_models = []

default_periodic_offset_start = -20
default_periodic_offset_end = 0
default_metagene_profile_length = 21

default_type_field = 'type'
default_position_field = 'position'
default_count_field = 'count'

default_iterations = 500
default_chains = 2
default_seed = 8675309

class NullDevice():
    def write(self, s):
        pass

def estimate_marginal_likelihoods(signal, periodic_models, nonperiodic_models, iterations, chains, seed):
        
    # construct the input for the models
    x_1 = signal[0::3]
    x_2 = signal[1::3]
    x_3 = signal[2::3]
    T = len(x_1)

    very_high_prior_location = max(signal)
    
    data = {
        "x_1": x_1,
        "x_2": x_2,
        "x_3": x_3,
        "T": T,
        "very_high_prior_location": very_high_prior_location
    }

    # get the likelihood for each of the models 
    bft_periodic = [
        pm.sampling(data=data, iter=iterations, chains=chains, n_jobs=1, seed=seed, refresh=-1)
            for pm in periodic_models]

    bft_nonperiodic = [
        nm.sampling(data=data, iter=iterations, chains=chains, n_jobs=1, seed=seed, refresh=-1)
            for nm in nonperiodic_models]
        
    return (bft_periodic, bft_nonperiodic)

def estimate_profile_bayes_factors(profile, args):
    length = profile['length'].iloc[0]

    # read in the relevant models
    periodic_models = [pickle.load(open(pm, 'rb')) for pm in args.periodic_models]
    nonperiodic_models = [pickle.load(open(npm, 'rb')) for npm in args.nonperiodic_models]
    
    # pull out the start offsets ("position" field) and counts
    mask_start = profile[args.type_field] == 'start'
    start_profile_df = profile.loc[mask_start]
    start_profile_df = start_profile_df.sort_values(args.position_field)

    start_positions = start_profile_df[args.position_field].values
    start_counts = start_profile_df[args.count_field].values
    
    # find the positions of the offsets of interest within the arrays
    begin_index = np.where(start_positions==args.periodic_offset_start)[0][0]
    stop_index = np.where(start_positions==args.periodic_offset_end)[0][0]
    
    # collect all of the results as a data frame
    ret = []
    
    for i in range(begin_index, stop_index+1):
        offset = start_positions[i]

        msg = "Length: {}, Offset: {}".format(length, offset)
        logger.debug(msg)

        # pull out the signal for this offset
        signal = start_counts[i:i+args.metagene_profile_length]
        (bft_periodic, bft_nonperiodic) = estimate_marginal_likelihoods(signal, 
            periodic_models, nonperiodic_models, 
            iterations=args.iterations,chains=args.chains,seed=args.seed)
    
        # extract the parameters of interest
        m_periodic_ex = [m.extract(pars=['lp__']) for m in bft_periodic]
        m_nonperiodic_ex = [m.extract(pars=['lp__']) for m in bft_nonperiodic]
        

        # now, choose the best model of each class,  based on mean likelihood
        m_periodic_means = [np.mean(m_ex['lp__']) for m_ex in m_periodic_ex]
        m_nonperiodic_means = [np.mean(m_ex['lp__']) for m_ex in m_nonperiodic_ex]

        max_periodic_mean = np.argmax(m_periodic_means)
        max_nonperiodic_mean = np.argmax(m_nonperiodic_means)

        # select the best sampling results
        m_periodic_ex = m_periodic_ex[max_periodic_mean]
        m_nonperiodic_ex = m_nonperiodic_ex[max_nonperiodic_mean]

        profile_sum = np.sum(signal)
        profile_peak = signal[0]

        v = {
            "offset": offset,
            "p_periodic_mean": np.mean(m_periodic_ex['lp__']),
            "p_periodic_var": np.var(m_periodic_ex['lp__']),
            "p_nonperiodic_mean": np.mean(m_nonperiodic_ex['lp__']),
            "p_nonperiodic_var": np.var(m_nonperiodic_ex['lp__']),
            'profile_sum': profile_sum,
            'profile_peak': profile_peak
        }

        v['bayes_factor_mean'] = v['p_periodic_mean'] - v['p_nonperiodic_mean']
        v['bayes_factor_var'] = v['p_periodic_var'] + v['p_periodic_var']

        ret.append(pd.Series(v))
        
    ret = pd.DataFrame(ret)
    ret['length'] = length

    return ret

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script estimates the Bayes factors for all metagene profiles in the "
        "given file. The script accepts as input multiple \"periodic\" and \"nonperiodic\" "
        "models. It uses the models of each type with the best mean to estimate the Bayes "
        "factor distributions.\n\nIt contains some hard-coded field names.")
        
    parser.add_argument('metagene_profiles', help="The (csv) file containing the metagene profiles")
    parser.add_argument('out', help="The output (csv.gz) file")

    parser.add_argument('--periodic-models', help="A list of pickled StanModel files which contain "
        "models that somehow represent periodic metagene profiles", nargs="+", 
        default=default_periodic_models)
    parser.add_argument('--nonperiodic-models', help="A list of pickled StanModel files which contain "
        "models that somehow represent nonperiodic metagene profiles", nargs="+", 
        default=default_nonperiodic_models)


    parser.add_argument('--periodic-offset-start', help="The position, relative to the translation "
        "initiation site, to begin calculating periodicity Bayes factors (inclusive)", type=int,
        default=default_periodic_offset_start)
    parser.add_argument('--periodic-offset-end', help="The position, relative to the translation "
        "initiation site, to stop calculating periodicity Bayes factors (inclusive)", type=int,
        default=default_periodic_offset_end)
    parser.add_argument('--metagene-profile-length', help="The length of the profile to use in the "
        "models. metagene_profile_length + periodic_offset_end must be consistent with the length "
        "of the extracted metagene profile. The length must be divisible by three.", type=int,
        default=default_metagene_profile_length)


    parser.add_argument('-s', '--seed', help="The random seeds to use for inference",
        type=int, default=default_seed)
    parser.add_argument('-c', '--chains', help="The number of MCMC chains to use", type=int,
        default=default_chains)
    parser.add_argument('-i', '--iterations', help="The number of MCMC iterations to use for "
        "each chain", type=int, default=default_iterations)
    
    parser.add_argument('-p', '--num-cpus', help="The number of CPUs to use. Each read "
        "length will be processed in its own thread (so that is the maximum number of CPUs "
        "that is useful).", type=int, default=default_num_cpus)

    parser.add_argument('--type-field', default=default_type_field)
    parser.add_argument('--count-field', default=default_count_field)
    parser.add_argument('--position-field', default=default_position_field)

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    # we will parallelize based on the lengths. So we need to know which lengths
    # are present in the metagene profiles file
    metagene_profiles = pd.read_csv(args.metagene_profiles)
    lengths = list(metagene_profiles['length'].unique())

    length_str = ','.join(str(int(l)) for l in lengths)
    msg = "Profiles will be created for lengths: {}".format(length_str)
    logger.info(msg)

    length_groups = metagene_profiles.groupby('length')

    # this does not seem to work
    # it is an attempt to redirect stderr to /dev/null so we don't
    # see the useless Stan output
    temp = sys.stderr
    f = open(os.devnull, 'w')
    sys.stderr = f

    all_profile_estimates_df = parallel.apply_parallel_groups(length_groups, args.num_cpus,
        estimate_profile_bayes_factors, args, progress_bar=True)
    sys.stderr = temp
    
    all_profile_estimates_df = pd.concat(all_profile_estimates_df)

    utils.write_df(all_profile_estimates_df, args.out, index=False)

if __name__ == '__main__':
    main()

