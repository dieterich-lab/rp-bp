#! /usr/bin/env python3

import logging
import os
import pandas as pd

import rpbp.ribo_utils.filenames as filenames

from rpbp.defaults import metagene_options, translation_options

logger = logging.getLogger(__name__)


class _return_key_dict(dict):
    def __missing__(self, key):
        return key


###
#   The following functions are helpful for parsing information out of the identifiers.
###


def get_transcript_id(orf_id, sep="_"):

    return orf_id.split(sep)[0]


def get_all_transcript_ids(orfs, sep="_", num_cpus=1, progress_bar=False):

    import pbiotools.misc.parallel as parallel

    transcript_ids = parallel.apply_parallel_iter(
        orfs["id"], num_cpus, get_transcript_id, sep, progress_bar=progress_bar
    )

    return transcript_ids


###
#   The following functions are all used for parsing replicates, etc., from the config file.
###


def get_riboseq_replicates(config):

    if "riboseq_biological_replicates" in config:
        if config["riboseq_biological_replicates"] is not None:
            msg = "Found 'riboseq_biological_replicates' key in config file"
            logger.info(msg)

            return config["riboseq_biological_replicates"]

    msg = (
        "Did not find 'riboseq_biological_replicates' key in config file. "
        "Using each 'riboseq_sample' as a single-condition replicate."
    )
    logger.info(msg)

    # create a dictionary mapping from the sample name to a single-element list
    ret = {name: [name] for name, sample in config["riboseq_samples"].items()}

    return ret


def get_riboseq_replicates_reverse_map(config):
    """Extract a mapping from sample to condition."""
    riboseq_replicates = get_riboseq_replicates(config)
    reverse_map = {v: k for k, l in riboseq_replicates.items() for v in l}

    ret_reverse_map = _return_key_dict()
    ret_reverse_map.update(reverse_map)

    return ret_reverse_map


def get_riboseq_cell_type_samples(config):
    if "riboseq_cell_type_samples" in config:
        if config["riboseq_cell_type_samples"] is not None:
            msg = "Found 'riboseq_cell_type_samples' key in config file"
            logger.info(msg)
            return config["riboseq_cell_type_samples"]

    msg = (
        "Did not find 'riboseq_cell_type_samples' key in config file. Using "
        "riboseq conditions (biological_replicate entries) as the cell types"
    )
    logger.info(msg)

    riboseq_replicates = get_riboseq_replicates(config)
    cell_type_samples = {x: [x] for x in riboseq_replicates}

    return cell_type_samples


def get_sample_name_map(config):
    """Extract the mapping from the 'riboseq_sample_name_map', or create
    a default one for all samples without an entry.
    """

    sample_name_map = _return_key_dict()

    if "riboseq_sample_name_map" in config:
        sample_name_map.update(config["riboseq_sample_name_map"])

    return sample_name_map


def get_condition_name_map(config):
    """Extract the mapping from the 'condition_name_map' and create a default
    one for all conditions without an entry.
    """

    condition_name_map = _return_key_dict()

    if "riboseq_condition_name_map" in config:
        condition_name_map.update(config["riboseq_condition_name_map"])

    return condition_name_map


###
#   The following functions are used for periodicity estimation, ORF filtering, etc.
###


def get_periodic_lengths_and_offsets(
    config,
    name,
    do_not_call=False,
    is_unique=True,
    default_params=None,
):

    """This function applies a set of filters to metagene profiles to select those
    which are "periodic" based on the read counts and Bayes factor estimates.

    First, the function checks if the configuration file sets the
    'use_fixed_lengths' flag is set. If so, then the specified lengths and
    offsets are returned.

    Otherwise, the function opens the appropriate file and extracts the filter
    values from the configuration file. In particular, it looks for the
    following keys:

    min_metagene_profile_count (float) : the minimum number of reads for a
        particular length in the filtered genome profile. Read lengths with
        fewer than this number of reads will not be used. default: 1000

    min_metagene_bf_mean (float) : if max_metagene_profile_bayes_factor_var
        is not None, then this is taken as a hard threshold on the estimated
        Bayes factor mean. If min_metagene_profile_bayes_factor_likelihood is
        given, then this is taken as the boundary value; that is, a profile is
        "periodic" if:

                [P(bf > min_metagene_bf_mean)] > min_metagene_bf_likelihood

        If both max_metagene_bf_var and min_metagene_bf_likelihood are None,
        then this is taken as a hard threshold on the mean for selecting
        periodic read lengths.

        If both max_metagene_bf_var and min_metagene_bf_likelihood are given,
        then both filters will be applied and the result will be the intersection.

    max_metagene_bf_var (float) : if given, then this is taken as a hard threshold
        on the estimated Bayes factor variance. default: None (i.e., this filter
        is not used)

    min_metagene_bf_likelihood (float) : if given, then this is taken a threshold
        on the likelihood of periodicity (see min_metagene_bf_mean description
        for more details). default: 0.5

    Parameters
    ----------
    config: dictionary
        the configuration information(see description)

    name: string
        the name of the dataset in question

    do_not_call: bool
        whether the metagene bf file should exist. If false, then dummy
        values are returned (and a warning message is printed).

    is_unique: bool
        whether only unique reads are used in the files

    default_params: default parameters, always passed when called from
        the Rp-Bp pipeline

    Returns
    -------
    lengths: list of strings
        all of the periodic read lengths

    offsets: list of strings
        the corresponding P-site offsets for the read lengths
    """
    import numpy as np
    import scipy.stats

    # check if we specified to just use a fixed offset and length
    if config.get("use_fixed_lengths", False):
        lengths = [str(l) for l in config["lengths"]]
        offsets = [str(o) for o in config["offsets"]]

        return (lengths, offsets)

    if default_params is None:
        default_params = metagene_options

    # filter out the lengths which do not satisfy the quality thresholds
    min_metagene_profile_count = config.get(
        "min_metagene_profile_count", default_params["min_metagene_profile_count"]
    )

    min_bf_mean = config.get(
        "min_metagene_bf_mean", default_params["min_metagene_bf_mean"]
    )

    max_bf_var = config.get(
        "max_metagene_bf_var", default_params["max_metagene_bf_var"]
    )

    min_bf_likelihood = config.get(
        "min_metagene_bf_likelihood", default_params["min_metagene_bf_likelihood"]
    )

    note_str = config.get("note", None)

    periodic_offsets = filenames.get_periodic_offsets(
        config["riboseq_data"],
        name,
        is_unique=is_unique,
        note=note_str,
    )

    if not os.path.exists(periodic_offsets):
        msg = (
            "The periodic offsets file does not exist. Please ensure the "
            "select-periodic-offsets script completed successfully or specify "
            'the "use_fixed_lengths", "lengths", and "offsets" values '
            "in the configuration file. '{}'".format(periodic_offsets)
        )

        if do_not_call:
            msg = msg + (
                '\nThe --do-not-call flag was given, so a "dummy" '
                "default length (29) and offset (12) will be used to check "
                "the remaining calls.\n"
            )

            logger.warning(msg)

            offsets = ["12"]
            lengths = ["29"]
            return (lengths, offsets)
        else:
            raise FileNotFoundError(msg)

    offsets_df = pd.read_csv(periodic_offsets)

    # we always use the count filter
    m_count = offsets_df["highest_peak_profile_sum"] > min_metagene_profile_count

    # which bf mean/variance filters do we use?
    m_bf_mean = True
    m_bf_var = True
    m_bf_likelihood = True

    if max_bf_var is not None:
        m_bf_mean = offsets_df["highest_peak_bf_mean"] > min_bf_mean
        m_bf_var = offsets_df["highest_peak_bf_var"] < max_bf_var

        msg = "Using the mean and variance filter. min_mean: {}, max_var: {}".format(
            min_bf_mean, max_bf_var
        )
        logger.debug(msg)

    if min_bf_likelihood is not None:
        # first, calculate the likelihood that the true BF is greater than m_bf_mean

        # the likelihood that BF>min_mean is 1-cdf(estimated_mean, estimated_var)

        # scipy parameterizes the normal using the std, so use sqrt(var)

        likelihood = 1 - scipy.stats.norm.cdf(
            min_bf_mean,
            offsets_df["highest_peak_bf_mean"],
            np.sqrt(offsets_df["highest_peak_bf_var"]),
        )

        nans = np.isnan(likelihood)
        num_nans = sum(nans)
        num_predictions = len(likelihood)

        max_likelihood = max(likelihood[~nans])

        # now filter
        m_bf_likelihood = likelihood > min_bf_likelihood

    if (max_bf_var is None) and (min_bf_likelihood is None):
        m_bf_mean = offsets_df["highest_peak_bf_mean"] > min_bf_mean

    filtered_periodic_offsets = offsets_df[
        m_count & m_bf_mean & m_bf_var & m_bf_likelihood
    ]

    offsets = filtered_periodic_offsets["highest_peak_offset"]
    lengths = filtered_periodic_offsets["length"]

    if len(lengths) == 0:
        msg = (
            "The periodic offsets file was found, but no periodic lengths "
            "were found. Please ensure the select-periodic-offsets script "
            'completed successfully or specify the "use_fixed_lengths", '
            '"lengths", and "offsets" values in the configuration file. '
            "'{}'".format(periodic_offsets)
        )

        if do_not_call:
            msg = msg + (
                '\nThe --do-not-call flag was given, so a "dummy" '
                "default length (29) and offset (12) will be used to check "
                "the remaining calls.\n"
            )

            logger.warning(msg)

            offsets = ["12"]
            lengths = ["29"]
            return (lengths, offsets)
        else:
            raise ValueError(msg)

    # offsets must be positive
    offsets = [str(-1 * int(o)) for o in offsets]
    lengths = [str(int(l)) for l in lengths]

    return (lengths, offsets)


def get_p_sites(bam_file, periodic_lengths, offsets):
    """Given a bam file of mapped riboseq reads, this function filters
    out the reads of non-periodic length, adjusts the start and end
    positions based on strand, and then shifts the remaining reads
    based on the length-specific offset.

    Args:
        bam_file (string) : the path to the mapped riboseq reads

        periodic_lengths (list-like) : a list of lengths to keep

        offsets (list-like) : the distance to shift each read of the
            respective length. the order here must match that in
            periodic_lengths

    Returns:
        pd.DataFrame : a data frame containing the transformed reads,
            sorted by chrom and start

    Imports:
        sys
        numpy
        pandas
        tqdm
        pysam
        bio_utils.bio
    """
    import sys
    import numpy as np
    import pandas as pd
    import tqdm

    import pysam
    import pbiotools.utils.bed_utils as bed_utils

    msg = "Reading BAM file"
    logger.info(msg)

    bam = pysam.AlignmentFile(bam_file)
    alignments = bam.fetch()
    num_alignments = bam.count()

    logger.info("Processing alignments")

    lengths = np.zeros(num_alignments, dtype=int)
    starts = np.zeros(num_alignments, dtype=int)
    ends = np.zeros(num_alignments, dtype=int)
    seqs = [""] * num_alignments
    strands = ["+"] * num_alignments
    fractions = np.zeros(num_alignments, dtype=float)

    al_iter = tqdm.tqdm(alignments, leave=True, file=sys.stdout, total=num_alignments)
    for i, a in enumerate(al_iter):
        starts[i] = a.reference_start
        ends[i] = a.reference_end
        lengths[i] = a.qlen
        seqs[i] = a.reference_name

        if a.is_reverse:
            strands[i] = "-"

    # The data frame will later be converted to BED6, so put the fields in the
    # correct order.
    map_df = pd.DataFrame()
    map_df["seqname"] = seqs
    map_df["start"] = starts
    map_df["end"] = ends
    map_df["id"] = "."
    map_df["score"] = "."
    map_df["strand"] = strands
    map_df["length"] = lengths

    msg = "Filtering reads by length"
    logger.info(msg)

    # now, filter based on lengths
    m_length = map_df["length"].isin(periodic_lengths)
    map_df = map_df[m_length]

    # now, we need to update the starts and ends based on the strand
    msg = "Updating coordinates based on offsets"
    logger.info(msg)

    # if the strand is positive, the end is start+1
    # if the strand is negative, the start is end-1
    m_positive = map_df["strand"] == "+"
    m_negative = map_df["strand"] == "-"

    # first, shift in the appropriate direction
    for i in range(len(periodic_lengths)):
        length = periodic_lengths[i]
        offset = offsets[i]

        m_length = map_df["length"] == length

        # adjust the start of forward strand
        map_df.loc[m_positive & m_length, "start"] = (
            map_df.loc[m_positive & m_length, "start"] + offset
        )

        # adjust the ends of negative strand
        map_df.loc[m_negative & m_length, "end"] = (
            map_df.loc[m_negative & m_length, "end"] - offset
        )

    # finally, we only care about the 5' end of the read, so discard everything else
    msg = "Discarding 3' end of reads"
    logger.info(msg)

    map_df.loc[m_positive, "end"] = map_df.loc[m_positive, "start"] + 1
    map_df.loc[m_negative, "start"] = map_df.loc[m_negative, "end"] - 1

    # now sort everything
    msg = "Sorting reads by coordinates"
    logger.info(msg)

    map_df = map_df.sort_values(["seqname", "start"])

    # and we only want the BED6 fields
    map_df = map_df[bed_utils.bed6_field_names]

    return map_df


def smooth_profile(
    profile,
    reweighting_iterations=translation_options["smoothing_reweighting_iterations"],
    fraction=translation_options["smoothing_fraction"],
):

    """This function smoothes the given ORF profile using the frame-specific
    approach. It assumes the profile is a dense numpy array and that any
    filtering due to differences of counts in reading frames, lengths, etc.,
    has already been performed.

    Please see the statsmodels.api.nonparametric.lowess documentation for
    more information about reweighting_iterations and fraction.

    Args:
        profile (np.array of numbers): an array containing the observed
            ORF profile. In principle, this could already be normalized.

        reweighting_iterations (int): the number of reweighting iterations

        fraction (float): the percentage of the signal to use for smooothing

    Returns:
        np.array: the smoothed profile

    Imports:
        statsmodels.api.nonparametric.lowess
    """
    import statsmodels.api as sm

    lowess = sm.nonparametric.lowess
    import numpy as np

    smoothed_profile = np.zeros_like(profile)

    # split the signal based on frame
    x_1 = profile[0::3]
    x_2 = profile[1::3]
    x_3 = profile[2::3]
    exog = np.arange(len(x_1))

    # x_1
    endog = x_1
    smoothed_x_1 = lowess(
        endog,
        exog,
        is_sorted=True,
        return_sorted=False,
        it=reweighting_iterations,
        frac=fraction,
    )

    # x_2
    endog = x_2
    smoothed_x_2 = lowess(
        endog,
        exog,
        is_sorted=True,
        return_sorted=False,
        it=reweighting_iterations,
        frac=fraction,
    )

    # x_3
    endog = x_3
    smoothed_x_3 = lowess(
        endog,
        exog,
        is_sorted=True,
        return_sorted=False,
        it=reweighting_iterations,
        frac=fraction,
    )

    smoothed_profile[0::3] = smoothed_x_1
    smoothed_profile[1::3] = smoothed_x_2
    smoothed_profile[2::3] = smoothed_x_3

    return smoothed_profile


def get_base_filter(
    bf,
    min_profile=translation_options["orf_min_profile_count"],
    min_length=translation_options["orf_min_length"],
):
    """This function extracts the ORFs from the BF dataframe which meet the
    minimum requirements to be considered for prediction. Namely, these
    requirements are:

        * The minimum sum across all reading frames exceeds the specified minimum
        * The length exceeds the specified minimum length
        * The number of reads in the first reading frame exceeds the number in
            either of the other two reading frames (though not necessarily the
            other two reading frames combined)

    Args:
        bf (pd.DataFrame): a data frame containing the relevant ORF information

        min_signal (int) : the minimum sum across all reading frames to consider
            an ORF as translated

        min_length (int) : the minimum length of ORF to consider

    Returns:
        boolean mask: a mask of the input data frame indicating all ORFs which
            meet the filtering criteria
    """

    if min_profile is None:
        m_profile = bf["profile_sum"] > 0
    else:
        m_profile = bf["profile_sum"] > min_profile

    m_length = bf["orf_len"] > min_length
    m_x1_gt_x2 = bf["x_1_sum"] > bf["x_2_sum"]
    m_x1_gt_x3 = bf["x_1_sum"] > bf["x_3_sum"]

    m_base = m_profile & m_length & m_x1_gt_x2 & m_x1_gt_x3
    return m_base


def get_bf_filter(
    bf,
    min_bf_mean=translation_options["min_bf_mean"],
    max_bf_var=translation_options["max_bf_var"],
    min_bf_likelihood=translation_options["min_bf_likelihood"],
):

    """This function applies filters to the Bayes factor estimates to find all
    ORFs which should be predicted as translated. This does not consider the
    length and profile sums, so this filter would need to be combined with
    the get_base_filter filter to find the true set of predicted ORFs.

    Args:
        bf (pd.DataFrame) : a data frame containing the relevant ORF information

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

    Returns:
        boolean mask: a mask of the input data frame indicating all ORFs which
            meet the filtering criteria

    Imports:
        numpy
        scipy.stats
    """
    import numpy as np
    import scipy.stats

    # which bf mean/variance filters do we use?
    m_bf_mean = True
    m_bf_var = True
    m_bf_likelihood = True

    if max_bf_var is not None:
        m_bf_mean = bf["bayes_factor_mean"] > min_bf_mean
        m_bf_var = bf["bayes_factor_var"] < max_bf_var
    if min_bf_likelihood is not None:
        # first, calculate the likelihood that the true BF is greater than m_bf_mean

        # the likelihood that BF>min_mean is 1-cdf(estimated_mean, estimated_var)

        # scipy parameterizes the normal using the std, so use sqrt(var)

        loc = bf["bayes_factor_mean"]
        scale = np.sqrt(bf["bayes_factor_var"])
        likelihood = 1 - scipy.stats.norm.cdf(min_bf_mean, loc, scale)

        nans = np.isnan(likelihood)
        num_nans = sum(nans)
        num_predictions = len(likelihood)

        msg = "Num nans: {}, num predictions: {}".format(num_nans, num_predictions)
        logger.debug(msg)

        if num_nans != num_predictions:
            max_likelihood = max(likelihood[~nans])
            msg = "Maximum likelihood: {}".format(max_likelihood)
            logger.debug(msg)

        # now filter
        m_bf_likelihood = likelihood > min_bf_likelihood

    if (max_bf_var is None) and (min_bf_likelihood is None):
        m_bf_mean = bf["bayes_factor_mean"] > min_bf_mean

    return m_bf_mean & m_bf_var & m_bf_likelihood


def get_predicted_orfs(
    bf,
    min_signal=translation_options["orf_min_profile_count"],
    min_length=translation_options["orf_min_length"],
    min_bf_mean=translation_options["min_bf_mean"],
    max_bf_var=translation_options["max_bf_var"],
    min_bf_likelihood=translation_options["min_bf_likelihood"],
    chisq_alpha=translation_options["chisq_alpha"],
    select_longest_by_stop=True,
    use_chi_square=False,
):
    """This function applies a set of filters to ORFs to select those which
    are predicted as "translated." This function selects translated ORFs
    based on the Bayes factor estimates or the chi-square p-values. ORFs
    must pass all of the relevant features to be selected as "translated."
    Optionally, among all ORFs which share a stop codon, only the longest
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

        select_longest_by_stop (bool): if True, then the selected ORFs will
            be merged based on stop codons: only the longest translated ORF
            at each stop codon will be returned. Otherwise, all ORFs will
            be returned.

        use_chi_square (bool): if True, then the selection is made based on
            the chi-square p-values only (Rp-chi), otherwise it is based on the Bayes
            factor estimates (Rp-Bp).

    Returns:
        all_orfs (pd.DataFrame) : all (longest) ORFs which meet the profile,
             length, frame filters

        predicted_orfs (pd.DataFrame) : all (longest) ORFs which meet the
            profile, length, frame Bayes factor (min_bf_mean, max_bf_var, min_bf_likelihood)
            or chisq_alpha filters

    Imports:
        bio_utils.bio
        numpy
        scipy.stats

    """
    import pbiotools.utils.bed_utils as bed_utils
    import scipy.stats

    msg = "Finding all ORFs with signal"
    logger.info(msg)

    m_base = get_base_filter(bf, min_signal, min_length)
    all_orfs = bf[m_base]

    # create the selected ORFs based on either Bayes factor or chisq_alpha
    if use_chi_square:
        M = len(all_orfs)
        # for the bonferroni correction, we only correct for the number of tests
        # we actually consider that is, we only correct for orfs which pass
        # the base filter
        corrected_significance_level = chisq_alpha / M

        msg = "Corrected significance level: {}".format(corrected_significance_level)
        logger.debug(msg)

        m_chisq_pval = all_orfs["chi_square_p"] < corrected_significance_level
        predicted_orfs = all_orfs[m_chisq_pval]
    else:
        m_bf = get_bf_filter(all_orfs, min_bf_mean, max_bf_var, min_bf_likelihood)
        predicted_orfs = all_orfs[m_bf]

    if select_longest_by_stop:
        all_orfs = bed_utils.get_longest_features_by_end(all_orfs)
        predicted_orfs = bed_utils.get_longest_features_by_end(predicted_orfs)

    return (all_orfs, predicted_orfs)
