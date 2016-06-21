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

    note_str = config.get('note', None)

    periodic_offsets = filenames.get_periodic_offsets(config['riboseq_data'], name, 
        is_unique=True, note=note_str)
    
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
    m_count = offsets_df['highest_peak_profile_sum'] > min_metagene_profile_count
    m_bf_mean = offsets_df['highest_peak_bf_mean'] > min_metagene_profile_bayes_factor_mean
    m_bf_var = offsets_df['highest_peak_bf_var'] < max_metagene_profile_bayes_factor_var

    filtered_periodic_offsets = offsets_df[m_count & m_bf_mean & m_bf_var]

    offsets = filtered_periodic_offsets['highest_peak_offset']
    lengths = filtered_periodic_offsets['length']

    # offsets must be positive
    offsets = [str(-1*int(o)) for o in offsets]
    lengths = [str(int(l)) for l in lengths]

    return (lengths, offsets)



default_min_signal = None
default_min_bf_mean = 5
default_max_bf_var = None
default_min_bf_likelihood = None
default_min_length = 0
default_chisq_alpha = 0.01

def get_predicted_orfs(bf, min_signal=default_min_signal, 
                            min_length=default_min_length,
                            min_bf_mean=default_min_bf_mean, 
                            max_bf_var=default_max_bf_var,
                            min_bf_likelihood=default_min_bf_likelihood,
                            chisq_alpha=default_chisq_alpha):
    """ This function applies a set of filters to ORFs to select those which
        are predicted as "translated." This function selects translated ORFs
        based on the Bayes factor estimates and the chi-square p-values. ORFs
        must pass all of the relevant features to be selected as "translated."
        Finally, among all ORFs which share a stop codon, only longest
        "translated" ORF is selected.

        Furthermore, for both BF and chi-square predictions, only ORFs which
        have more reads in the first reading frame than either of the other two
        will be selected as translated. (This is called the 'frame filter'
        below.)

        Args:
            bf (pd.DataFrame) : a data frame containing the relevant ORF information

            min_signal (int) : the minimum sum across all reading frames to consider
                an ORF as translated
            
            min_length (int) : the minimum length of ORF to consider

            min_bf_mean (float) : if max_bf_var is not None, then this is taken
                as a hard threshold on the estimated Bayes factor mean. If
                min_bf_likelihood is given, then this is taken as the boundary
                value; that is, an ORF is "translated" if:

                    [P(bf > min_bf_mean)] > min_bf_likelihood

                If both max_bf_var and min_bf_likelihood are None, then this is
                taken as a hard threshold on the mean for selecting translated ORFs.

                If both max_bf_var and min_bf_likelihood are given, then both
                filters will be applied and the result will be the intersection.

            max_bf_var (float) : if given, then this is taken as a hard threshold
                on the estimated Bayes factor variance

            min_bf_likelihood (float) : if given, then this is taken a threshold
                on the likelihood of translation (see min_bf_mean description
                for more details)

            chisq_alpha (float) : the significance value for selecting translated
                ORFs according to the chi-square test. This value is 
                Bonferroni-corrected based on the number of ORFs which meet the
                length, profile and frame filters.

        Returns:
            longest_orfs (pd.DataFrame) : all longest ORFs which meet the profile,
                 length, frame filters

            bf_longest_orfs (pd.DataFrame) : all longest ORFs which meet the
                profile, length, frame (min_bf_mean, max_bf_var, min_bf_likelihood) filters

            chisq_longest_orfs (pd.DataFrame) : all longest ORFs which meet the
                profile, length, frame, chisq_alpha filters

        Imports:
            misc.bio
            numpy
            scipy.stats

    """

    import misc.bio as bio
    import numpy as np
    import scipy.stats

    msg = "Finding all longest ORFs with signal"
    logging.info(msg)

    if min_signal is None:
        m_profile = bf['profile_sum'] > 0
    else:
        m_profile = bf['profile_sum'] > min_signal

    m_length = bf['orf_len'] > min_length
    m_x1_gt_x2 = bf['x_1_sum'] > bf['x_2_sum']
    m_x1_gt_x3 = bf['x_1_sum'] > bf['x_3_sum']

    m_base = m_profile & m_length & m_x1_gt_x2 & m_x1_gt_x3

    longest_orfs = bio.get_longest_features_by_end(bf[m_base])
    
    # create the selected ORFs

    # which bf mean/variance filters do we use? 
    m_bf_mean = True
    m_bf_var = True
    m_bf_likelihood = True

    if max_bf_var is not None:
        m_bf_mean = bf['bayes_factor_mean'] > min_bf_mean
        m_bf_var = bf['bayes_factor_var'] < max_bf_var
    if min_bf_likelihood is not None:
        # first, calculate the likelihood that the true BF is greater than m_bf_mean

        # the likelihood that BF>min_mean is 1-cdf(estimated_mean, estimated_var)

        # scipy parameterizes the normal using the std, so use sqrt(var)

        likelihood = 1-scipy.stats.norm.cdf(min_bf_mean, bf['bayes_factor_mean'], np.sqrt(bf['bayes_factor_var']))

        nans = np.isnan(likelihood)
        num_nans = sum(nans)
        num_predictions = len(likelihood)

        msg = "Num nans: {}, num predictions: {}".format(num_nans, num_predictions)
        logging.debug(msg)

        max_likelihood = max(likelihood[~nans])
        msg = "Maximum likelihood: {}".format(max_likelihood)
        logging.debug(msg)

        # now filter
        m_bf_likelihood = likelihood > min_bf_likelihood

    if (max_bf_var is None) and (min_bf_likelihood is None):
        m_bf_mean = bf['bayes_factor_mean'] > min_bf_mean

    # apply all the filters
    m_bf_predicted = m_base & m_bf_mean & m_bf_var & m_bf_likelihood

    bf_longest_predicted_orfs = bio.get_longest_features_by_end(bf[m_bf_predicted])

    M = len(longest_orfs)
    # for the bonferroni correction, we only correct for the number of tests we actually consider
    # that is, we only correct for orfs which pass the base filter
    corrected_significance_level = chisq_alpha / M

    msg = "Corrected significance level: {}".format(corrected_significance_level)
    logging.debug(msg)
    
    m_chisq_pval = bf['chi_square_p'] < corrected_significance_level
    m_chisq_predicted = m_base & m_chisq_pval

    chisq_longest_predicted_orfs = bio.get_longest_features_by_end(bf[m_chisq_predicted])
    
    return (longest_orfs, bf_longest_predicted_orfs, chisq_longest_predicted_orfs)

