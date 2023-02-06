#! /usr/bin/env python3

import sys
import logging
import argparse
import ctypes
import multiprocessing
import scipy.io
import scipy.stats
import scipy.sparse
import tempfile
import pathlib

import numpy as np
import pandas as pd

import pbiotools.utils.bed_utils as bed_utils
import pbiotools.misc.logging_utils as logging_utils
import pbiotools.misc.parallel as parallel
import pbiotools.misc.slurm as slurm
import pbiotools.misc.utils as utils

import rpbp.ribo_utils.utils as ribo_utils
import rpbp.ribo_utils.compile_rpbp_models as compile_rpbp_models

from cmdstanpy import CmdStanModel
from rpbp.defaults import default_num_groups, translation_options

logger = logging.getLogger(__name__)
cmdstanpy_logger = logging.getLogger("cmdstanpy")

# we will use global variables to share the (read-only) scipy.sparse.csr_matrix
# across the child processes.

# see: http://stackoverflow.com/questions/1675766/
# how-to-combine-pool-map-with-array-shared-memory-in-python-multiprocessing
profiles_data = 0
profiles_indices = 0
profiles_indptr = 0
profiles_shape = 0

translated_models = 0
untranslated_models = 0
args = 0


def get_bayes_factor(profile, translated_models, untranslated_models, args):
    """This function calculates the Bayes' factor for a single ORF profile.

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
    """
    profile_sum = sum(profile)

    # split the signal based on frame
    x_1 = profile[0::3]
    x_2 = profile[1::3]
    x_3 = profile[2::3]
    T = len(x_1)

    x_1_sum = sum(x_1)
    x_2_sum = sum(x_2)
    x_3_sum = sum(x_3)

    ret = {
        "p_translated_mean": float("-inf"),
        "p_translated_var": float("-inf"),
        "p_background_mean": float("-inf"),
        "p_background_var": float("-inf"),
        "bayes_factor_mean": float("-inf"),
        "bayes_factor_var": float("-inf"),
        "chi_square_p": float("-inf"),
        "x_1_sum": x_1_sum,
        "x_2_sum": x_2_sum,
        "x_3_sum": x_3_sum,
        "profile_sum": profile_sum,
    }
    ret = pd.Series(ret)

    # check if something odd happens with the length
    # this should already be checked before calling the function.
    if (T != len(x_2)) or (T != len(x_3)):
        return ret

    # and make sure we have more reads in x_1 than each of the others
    if (x_1_sum < x_2_sum) or (x_1_sum < x_3_sum):
        return ret

    # chi-square values
    f_obs = [x_1_sum, x_2_sum, x_3_sum]
    chisq, chi_square_p = scipy.stats.chisquare(f_obs)
    ret["chi_square_p"] = chi_square_p

    # check if we only wanted the chi square value
    if args.chi_square_only:
        return ret

    # now, smooth the signals
    smoothed_profile = ribo_utils.smooth_profile(
        profile,
        reweighting_iterations=args.reweighting_iterations,
        fraction=args.fraction,
    )

    # split the signal based on frame
    x_1 = smoothed_profile[0::3]
    x_2 = smoothed_profile[1::3]
    x_3 = smoothed_profile[2::3]
    nonzero_x_1 = np.count_nonzero(x_1)

    # construct the input for Stan
    data = {"x_1": x_1, "x_2": x_2, "x_3": x_3, "T": T, "nonzero_x_1": nonzero_x_1}

    iter_warmup = int(args.iterations // 2)

    with tempfile.TemporaryDirectory() as tmpdir:
        m_translated = [
            tm.sample(
                data=data,
                iter_warmup=iter_warmup,
                iter_sampling=iter_warmup,
                chains=args.chains,
                parallel_chains=1,
                seed=args.seed,
                show_progress=False,
                show_console=False,
                output_dir=tmpdir,
            )
            for tm in translated_models
        ]

        m_background = [
            bm.sample(
                data=data,
                iter_warmup=iter_warmup,
                iter_sampling=iter_warmup,
                chains=args.chains,
                parallel_chains=1,
                seed=args.seed,
                show_progress=False,
                show_console=False,
                output_dir=tmpdir,
            )
            for bm in untranslated_models
        ]

        # extract the parameters of interest
        m_translated_ex = [
            m.draws_pd()[["lp__", "signal_location", "signal_scale"]]
            for m in m_translated
        ]
        m_background_ex = [
            m.draws_pd()[["lp__", "background_location", "background_scale"]]
            for m in m_background
        ]

    # now, choose the best model of each class,  based on mean likelihood
    m_translated_means = [np.mean(m_ex["lp__"]) for m_ex in m_translated_ex]
    m_background_means = [np.mean(m_ex["lp__"]) for m_ex in m_background_ex]

    max_translated_mean = np.argmax(m_translated_means)
    max_background_mean = np.argmax(m_background_means)

    # select the best sampling results
    m_translated_ex = m_translated_ex[max_translated_mean]
    m_background_ex = m_background_ex[max_background_mean]

    # extract the relevant means and variances
    ret["p_translated_mean"] = np.mean(m_translated_ex["lp__"])
    ret["p_translated_var"] = np.var(m_translated_ex["lp__"])

    ret["p_background_mean"] = np.mean(m_background_ex["lp__"])
    ret["p_background_var"] = np.var(m_background_ex["lp__"])

    # the (log of) the Bayes factor is the difference between two normals:
    # (the best translated model) - (the best background model)
    #
    # thus, it is also a normal whose mean is the difference of the two means
    # and whose variance is the sum of the two variances
    ret["bayes_factor_mean"] = ret["p_translated_mean"] - ret["p_background_mean"]
    ret["bayes_factor_var"] = ret["p_translated_var"] + ret["p_background_var"]

    return ret


def get_all_bayes_factors(orfs):

    """This function calculates the Bayes' factor term for each region in regions. See the
    description of the script for the Bayes' factor calculations.

    Args:
        orfs (pd.DataFrame) : a set of orfs. The columns must include:
            orf_num
            exon_lengths

        args (namespace) : a namespace containing the models and profiles filenames

    Returns:
        pandas.Series: the Bayes' factors (and other estimated quantities) for each region
    """

    profiles = scipy.sparse.csr_matrix(
        (profiles_data, profiles_indices, profiles_indptr),
        shape=profiles_shape,
        copy=False,
    )

    logger.debug("Applying on regions")
    bfs = []
    for idx, row in orfs.iterrows():
        orf_num = row[args.orf_num_field]
        orf_len = row["orf_len"]

        # sometimes the orf_len is off...
        if orf_len % 3 != 0:
            msg = "Found an ORF whose length was not 0 mod 3. Skipping. orf_id: {}".format(
                row["id"]
            )
            logger.warning(msg)
            continue

        profile = utils.to_dense(profiles, orf_num, float, length=orf_len)

        row_bf = get_bayes_factor(profile, translated_models, untranslated_models, args)
        row = pd.concat([row, row_bf])

        bfs.append(row)

    bfs = pd.DataFrame(bfs)
    return bfs


def main():
    global profiles_data, profiles_indices, profiles_indptr, profiles_shape
    global translated_models, untranslated_models
    global args

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""This script uses Hamiltonian MCMC with Stan
        to estimate translation parameters for a set of regions (presumably ORFs). Roughly, it takes
        as input: (1) a set of regions (ORFs) and their corresponding profiles
                  (2) a "translated" model which gives the probability that a region is translated
                  (3) an "untranslated" model which gives the probability that a region is not translated.
        The script first smoothes the profiles using LOWESS. It then calculates both the Bayes' factor
        (using the smoothed profile) and chi2 value (using the raw counts) for each ORF.""",
    )

    parser.add_argument("profiles", help="The ORF profiles (counts) (mtx)")

    parser.add_argument(
        "regions", help="The regions (ORFs) for which predictions will be made (BED12+)"
    )

    parser.add_argument("out", help="The output file for the Bayes' factors (BED12+)")

    parser.add_argument(
        "--chi-square-only",
        help="""If this flag is present, then only the chi
        square test will be performed for each ORF. This can also be a way to get the counts within
        each of the ORFs.""",
        action="store_true",
    )

    parser.add_argument(
        "--translated-models", help="The models to use as H_t (pkl)", nargs="+"
    )

    parser.add_argument(
        "--untranslated-models", help="The models to use as H_u (pkl)", nargs="+"
    )

    # filtering options
    parser.add_argument(
        "--min-length",
        help="ORFs with length less than this value will not be processed",
        type=int,
        default=translation_options["orf_min_length_pre"],
    )

    parser.add_argument(
        "--max-length",
        help="ORFs with length greater than this value will not be processed",
        type=int,
        default=translation_options["orf_max_length_pre"],
    )

    parser.add_argument(
        "--min-profile",
        help="""ORFs with profile sum (i.e., number of reads) less than this
        value will not be processed.""",
        type=float,
        default=translation_options["orf_min_profile_count_pre"],
    )

    # smoothing options
    parser.add_argument(
        "--fraction",
        help="The fraction of signal to use in LOWESS",
        type=float,
        default=translation_options["smoothing_fraction"],
    )

    parser.add_argument(
        "--reweighting-iterations",
        help="The number of reweighting "
        "iterations to use in LOWESS. "
        "Please see the statsmodels documentation for a "
        "detailed description of this parameter.",
        type=int,
        default=translation_options["smoothing_reweighting_iterations"],
    )

    # MCMC options
    parser.add_argument(
        "-s",
        "--seed",
        help="The random seeds to use for inference",
        type=int,
        default=translation_options["seed"],
    )
    parser.add_argument(
        "-c",
        "--chains",
        help="The number of MCMC chains to use",
        type=int,
        default=translation_options["chains"],
    )
    parser.add_argument(
        "-i",
        "--iterations",
        help="The number of MCMC iterations to use for each chain",
        type=int,
        default=translation_options["translation_iterations"],
    )

    # behavior options
    # [--num-orfs] is not used in the the Rp-Bp pipeline
    parser.add_argument(
        "--num-orfs",
        help="If n>0, then only this many ORFs will be processed",
        type=int,
        default=0,
    )

    parser.add_argument("--orf-num-field", default="orf_num")

    parser.add_argument(
        "--do-not-compress",
        help="Unless otherwise specified, the output will " "be written in GZip format",
        action="store_true",
    )

    parser.add_argument(
        "-g",
        "--num-groups",
        help="The number of groups into which to split "
        "the ORFs. More groups means the progress bar is "
        "updated more frequently but incurs more overhead "
        "because of the parallel calls.",
        type=int,
        default=default_num_groups,
    )

    slurm.add_sbatch_options(parser)
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    if args.use_slurm:
        cmd = " ".join(sys.argv)
        slurm.check_sbatch(cmd, args=args)
        return

    # logging
    cmdstanpy_logger.disabled = True
    if args.enable_ext_logging:
        cmdstanpy_logger.disabled = False

    # read in the regions and apply the filters
    msg = "Reading and filtering ORFs"
    logger.info(msg)
    regions = bed_utils.read_bed(args.regions, low_memory=False)

    # by default, keep everything
    m_filters = np.array([True] * len(regions))

    # min length
    if args.min_length > 0:
        m_min_length = regions["orf_len"] >= args.min_length
        m_filters = m_min_length & m_filters

    # max length
    if args.max_length > 0:
        m_max_length = regions["orf_len"] <= args.max_length
        m_filters = m_max_length & m_filters

    # min profile
    profiles = scipy.io.mmread(args.profiles).tocsr()
    profiles_sums = profiles.sum(axis=1)
    good_orf_nums = np.where(profiles_sums >= args.min_profile)
    good_orf_nums = set(good_orf_nums[0])
    m_profile = regions["orf_num"].isin(good_orf_nums)
    m_filters = m_profile & m_filters

    regions = regions[m_filters]

    if args.num_orfs > 0:
        regions = regions.head(args.num_orfs)

    regions = regions.reset_index(drop=True)

    msg = "Number of regions after filtering: {}".format(len(regions))
    logger.info(msg)

    compile_rpbp_models.compile()

    # read in the relevant pre-compiled models
    translated_models = [
        CmdStanModel(exe_file=pathlib.Path(tm).with_suffix(""))
        for tm in args.translated_models
    ]
    untranslated_models = [
        CmdStanModel(exe_file=pathlib.Path(bm).with_suffix(""))
        for bm in args.untranslated_models
    ]

    profiles_data = multiprocessing.RawArray(ctypes.c_double, profiles.data.flat)
    profiles_indices = multiprocessing.RawArray(ctypes.c_int, profiles.indices)
    profiles_indptr = multiprocessing.RawArray(ctypes.c_int, profiles.indptr)
    profiles_shape = multiprocessing.RawArray(ctypes.c_int, profiles.shape)

    bfs_l = parallel.apply_parallel_split(
        regions,
        args.num_cpus,
        get_all_bayes_factors,
        num_groups=args.num_groups,
        progress_bar=True,
        backend="multiprocessing",
    )

    bfs = pd.concat(bfs_l)

    # write the results as a bed12+ file
    bed_utils.write_bed(bfs, args.out)


if __name__ == "__main__":
    main()
