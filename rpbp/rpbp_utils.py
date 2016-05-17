import logging
import os
import pandas as pd

import rpbp.filenames as filenames

default_min_metagene_profile_count = 1000
default_min_metagene_profile_bayes_factor_mean = 5
default_max_metagene_profile_bayes_factor_var = 5

def get_periodic_lengths_and_offsets(config, name, do_not_call=False):
    # check if we specified to just use a fixed offset and length
    if 'use_fixed_lengths' in config:
        lengths = config['lengths']
        offsets = config['offsets']

        return (lengths, offsets)

    # filter out the lengths which do not satisfy the quality thresholds
    min_metagene_profile_count = config.get(
        "min_metagene_profile_count", default_min_metagene_profile_count)

    min_metagene_profile_bayes_factor_mean = config.get(
        "min_metagene_profile_bayes_factor_mean", default_min_metagene_profile_bayes_factor_mean)

    max_metagene_profile_bayes_factor_var = config.get(
        "max_metagene_profile_bayes_factor_var", default_max_metagene_profile_bayes_factor_var)

    periodic_offsets = filenames.get_periodic_offsets(config['riboseq_data'], name, 
        is_unique=True)
    
    if not os.path.exists(periodic_offsets):
        msg = ("The periodic offsets file does not exist. Please ensure the select-periodic-offsets "
            "script completed successfully or specify the \"use_fixed_lengths\", \"lengths\", and "
            "\"offsets\" values in the configuration file. '{}'".format(periodic_offsets))

        if do_not_call:
            msg = msg +  ("\nThe --do-not-call flag was given, so \"dummy\" default lengths will be "
                "used to check the remaining calls.\n")

            logging.warning(msg)

            offsets = ["12"]
            lengths = ["29"]
            return (lengths, offsets)
        else:
            raise FileNotFoundError(msg)
    
    offsets_df = pd.read_csv(periodic_offsets)
    m_count = offsets_df['highest_peak_peak'] > min_metagene_profile_count
    m_bf_mean = offsets_df['highest_peak_bf_mean'] > min_metagene_profile_bayes_factor_mean
    m_bf_var = offsets_df['highest_peak_bf_var'] < max_metagene_profile_bayes_factor_var

    filtered_periodic_offsets = offsets_df[m_count & m_bf_mean & m_bf_var]

    offsets = filtered_periodic_offsets['highest_peak_offset']
    lengths = filtered_periodic_offsets['length']

    # offsets must be positive
    offsets = [str(-1*int(o)) for o in offsets]
    lengths = [str(int(l)) for l in lengths]

    return (lengths, offsets)



default_min_signal = 10
default_min_bf_mean = 5
default_max_bf_var = 2
default_min_length = 200
default_chisq_alpha = 0.01

def get_predicted_orfs(bf, min_signal=default_min_signal, min_bf_mean=default_min_bf_mean, 
                       max_bf_var=default_max_bf_var, min_length=default_min_length,
                       chisq_alpha=default_chisq_alpha):

    import misc.bio as bio

    msg = "Finding all longest ORFs with signal"
    logging.info(msg)

    mask_profile = bf['profile_sum'] > min_signal
    longest_orfs = bio.get_longest_features_by_end(bf[mask_profile])
    
    # create the selected ORFs
    m_profile_sum = bf['profile_sum'] > min_signal
    m_bf_mean = bf['bayes_factor_mean'] > min_bf_mean
    m_bf_var = bf['bayes_factor_var'] < max_bf_var
    m_x1_gt_x2 = bf['x_1_sum'] > bf['x_2_sum']
    m_x1_gt_x3 = bf['x_1_sum'] > bf['x_3_sum']
    
    m_length = bf['orf_len'] > min_length

    m_bf_predicted = m_x1_gt_x2 & m_x1_gt_x3 & m_profile_sum & m_bf_mean & m_bf_var & m_length

    bf_longest_predicted_orfs = bio.get_longest_features_by_end(bf[m_bf_predicted])

    M = len(longest_orfs)
    # for the bonferroni correction, we only correct for the number of tests we actually consider
    # that is, we only correct for orfs which pass the minimum profile filter
    corrected_significance_level = chisq_alpha / sum(m_profile_sum)

    msg = "Corrected significance level: {}".format(corrected_significance_level)
    logging.debug(msg)
    
    m_chisq_pval = bf['chi_square_p'] < corrected_significance_level
    m_chisq_predicted = m_x1_gt_x2 & m_x1_gt_x3 & m_profile_sum & m_chisq_pval & m_length

    chisq_longest_predicted_orfs = bio.get_longest_features_by_end(bf[m_chisq_predicted])
    
    return (longest_orfs, bf_longest_predicted_orfs, chisq_longest_predicted_orfs)

