#! /usr/bin/env python3

import argparse
import pandas as pd

import misc.utils as utils

def get_most_periodic_offset(profile_df):
    
    #mask_largest_count = profile_df['profile_sum'] == profile_df['profile_sum'].max()
    #largest_count_row = profile_df[mask_largest_count].iloc[0]
    m_highest_peak = profile_df['profile_peak'] == profile_df['profile_peak'].max()
    highest_peak_row = profile_df[m_highest_peak].iloc[0]

    length = highest_peak_row['length']

    highest_peak_peak = highest_peak_row['profile_peak']
    highest_peak_profile_sum = highest_peak_row['profile_sum']
    highest_peak_offset = int(highest_peak_row['offset'])
    highest_peak_bf_mean = highest_peak_row['bayes_factor_mean']
    highest_peak_bf_var = highest_peak_row['bayes_factor_var']

    # and now the offset with the best bayes factor mean
    mask_largest_bf = profile_df['bayes_factor_mean'] == profile_df['bayes_factor_mean'].max()
    largest_bf_row = profile_df[mask_largest_bf].iloc[0]
    
    largest_bf_peak = largest_bf_row['profile_peak']
    largest_bf_profile_sum = largest_bf_row['profile_sum']
    largest_bf_offset = int(largest_bf_row['offset'])
    largest_bf_mean = largest_bf_row['bayes_factor_mean']
    largest_bf_var = largest_bf_row['bayes_factor_var']

    # and construct the output
    ret = {
        'length': length,
        'highest_peak_peak': highest_peak_peak,
        'highest_peak_profile_sum': highest_peak_profile_sum,
        'highest_peak_offset': highest_peak_offset,
        'highest_peak_bf_mean': highest_peak_bf_mean,
        'highest_peak_bf_var': highest_peak_bf_var,
        'largest_bf_peak': largest_bf_peak,
        'largest_bf_profile_sum': largest_bf_profile_sum,
        'largest_bf_offset': largest_bf_offset,
        'largest_bf_mean': largest_bf_mean,
        'largest_bf_var': largest_bf_var
    }

    return pd.Series(ret)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script combines the Bayes' factors for the metagene profiles' "
        "periodicity and signal information to select the offsets for beginning the riboseq "
        "signal.\n\nThis script contains some hard-coded field names.")
    parser.add_argument('profile_bayes_factors', help="The (csv) file containing the Bayes' "
        "factors for periodicity at each length")
    parser.add_argument('out', help="The output (csv.gz) file")
    args = parser.parse_args()

    # read in the Bayes' factors and profiles
    bayes_factors_df = pd.read_csv(args.profile_bayes_factors)
    length_groups = bayes_factors_df.groupby('length')
    periodic_offsets = length_groups.apply(get_most_periodic_offset)
    utils.write_df(periodic_offsets, args.out, index=False)

if __name__ == '__main__':
    main()
