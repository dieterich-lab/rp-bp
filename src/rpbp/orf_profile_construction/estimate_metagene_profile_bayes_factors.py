#! /usr/bin/env python3

import argparse
import logging
import pathlib

import numpy as np
import pandas as pd

from cmdstanpy import CmdStanModel

import rpbp.ribo_utils.compile_rpbp_models as compile_rpbp_models

import pbiotools.misc.logging_utils as logging_utils
import pbiotools.misc.parallel as parallel
import pbiotools.misc.utils as utils
import pbiotools.misc.pandas_utils as pandas_utils

logger = logging.getLogger(__name__)
cmdstanpy_logger = logging.getLogger("cmdstanpy")

# When called from the Rp-Bp pipeline (create-orf-profiles), default
# options (or else specified via the configuration file) are always
# passed as arguments.

default_periodic_offset_start = -20
default_periodic_offset_end = 0
default_metagene_profile_length = 21

default_iterations = 500
default_chains = 2
default_seed = 8675309

# Not passed as arguments, unlikely to be required

default_type_field = "type"
default_position_field = "position"
default_count_field = "count"


def estimate_marginal_likelihoods(
    signal, periodic_models, nonperiodic_models, iterations, chains, seed
):

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
        "very_high_prior_location": very_high_prior_location,
    }

    iter_warmup = int(iterations // 2)

    # get the likelihood for each of the models
    bft_periodic = [
        pm.sample(
            data=data,
            iter_warmup=iter_warmup,
            iter_sampling=iter_warmup,
            chains=chains,
            parallel_chains=1,
            seed=seed,
            show_progress=False,
            show_console=False,
        )
        for pm in periodic_models
    ]

    bft_nonperiodic = [
        nm.sample(
            data=data,
            iter_warmup=iter_warmup,
            iter_sampling=iter_warmup,
            chains=chains,
            parallel_chains=1,
            seed=seed,
            show_progress=False,
            show_console=False,
        )
        for nm in nonperiodic_models
    ]

    return (bft_periodic, bft_nonperiodic)


def estimate_profile_bayes_factors(profile, args):

    # logging
    cmdstanpy_logger.disabled = True
    if args.enable_ext_logging:
        cmdstanpy_logger.disabled = False

    length = profile["length"].iloc[0]

    # read in the relevant pre-compiled models
    periodic_models = [
        CmdStanModel(exe_file=pathlib.Path(pm).with_suffix(""))
        for pm in args.periodic_models
    ]
    nonperiodic_models = [
        CmdStanModel(exe_file=pathlib.Path(npm).with_suffix(""))
        for npm in args.nonperiodic_models
    ]

    # pull out the start offsets ("position" field) and counts
    mask_start = profile[args.type_field] == "start"
    start_profile_df = profile.loc[mask_start]
    start_profile_df = start_profile_df.sort_values(args.position_field)

    start_positions = start_profile_df[args.position_field].values
    start_counts = start_profile_df[args.count_field].values

    # find the positions of the offsets of interest within the arrays
    begin_start_pos = np.where(start_positions == args.periodic_offset_start)[0]
    if len(begin_start_pos) == 0:
        msg = "Did not find any start offsets for length: {}".format(length)
        logging.warning(msg)
        return None
    begin_index = begin_start_pos[0]

    stop_start_pos = np.where(start_positions == args.periodic_offset_end)[0]
    if len(stop_start_pos) == 0:
        msg = "Did not find any stop offsets for length: {}".format(length)
        logging.warning(msg)
        return None
    stop_index = stop_start_pos[0]

    # collect all of the results as a data frame
    ret = []

    for i in range(begin_index, stop_index + 1):
        offset = start_positions[i]

        msg = "Length: {}, Offset: {}".format(length, offset)
        logger.debug(msg)

        # pull out the signal for this offset
        signal = start_counts[i : i + args.metagene_profile_length]
        (bft_periodic, bft_nonperiodic) = estimate_marginal_likelihoods(
            signal,
            periodic_models,
            nonperiodic_models,
            iterations=args.iterations,
            chains=args.chains,
            seed=args.seed,
        )

        # extract the parameters of interest
        m_periodic_ex = [m.draws_pd()["lp__"].values for m in bft_periodic]
        m_nonperiodic_ex = [m.draws_pd()["lp__"].values for m in bft_nonperiodic]

        # now, choose the best model of each class,  based on mean likelihood
        m_periodic_means = [np.mean(m_ex) for m_ex in m_periodic_ex]
        m_nonperiodic_means = [np.mean(m_ex) for m_ex in m_nonperiodic_ex]

        max_periodic_mean = np.argmax(m_periodic_means)
        max_nonperiodic_mean = np.argmax(m_nonperiodic_means)

        # select the best sampling results
        m_periodic_ex = m_periodic_ex[max_periodic_mean]
        m_nonperiodic_ex = m_nonperiodic_ex[max_nonperiodic_mean]

        profile_sum = np.sum(signal)
        profile_peak = signal[0]

        v = {
            "offset": offset,
            "p_periodic_mean": np.mean(m_periodic_ex),
            "p_periodic_var": np.var(m_periodic_ex),
            "p_nonperiodic_mean": np.mean(m_nonperiodic_ex),
            "p_nonperiodic_var": np.var(m_nonperiodic_ex),
            "profile_sum": profile_sum,
            "profile_peak": profile_peak,
        }

        v["bayes_factor_mean"] = v["p_periodic_mean"] - v["p_nonperiodic_mean"]
        v["bayes_factor_var"] = v["p_periodic_var"] + v["p_nonperiodic_var"]

        ret.append(pd.Series(v))

    ret = pd.DataFrame(ret)
    ret["length"] = length

    return ret


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""This script estimates the Bayes factors
        for all metagene profiles in the given file. The script accepts as input multiple
        'periodic' and 'nonperiodic' models. It uses the models of each type with the
        best mean to estimate the Bayes factor distributions. It contains some
        hard-coded field names.""",
    )

    parser.add_argument(
        "metagene_profiles", help="The (csv) file containing the metagene profiles"
    )

    parser.add_argument("out", help="The output (csv.gz) file")

    parser.add_argument(
        "--periodic-models",
        help="""A list of pickled StanModel files which contain
        models that somehow represent periodic metagene profiles.""",
        nargs="+",
        default=[],
    )

    parser.add_argument(
        "--nonperiodic-models",
        help=""""A list of pickled StanModel files which contain
        models that somehow represent nonperiodic metagene profiles.""",
        nargs="+",
        default=[],
    )

    parser.add_argument(
        "--periodic-offset-start",
        help="""The position, relative to the translation
        initiation site, to begin calculating periodicity Bayes factors (inclusive).""",
        type=int,
        default=default_periodic_offset_start,
    )

    parser.add_argument(
        "--periodic-offset-end",
        help="""The position, relative to the translation
        initiation site, to stop calculating periodicity Bayes factors (inclusive)""",
        type=int,
        default=default_periodic_offset_end,
    )

    parser.add_argument(
        "--metagene-profile-length",
        help="""The length of the profile to use in the
        models. metagene_profile_length + periodic_offset_end must be consistent with the length
        of the extracted metagene profile. The length must be divisible by three.""",
        type=int,
        default=default_metagene_profile_length,
    )

    parser.add_argument(
        "-s",
        "--seed",
        help="The random seeds to use for inference",
        type=int,
        default=default_seed,
    )

    parser.add_argument(
        "-c",
        "--chains",
        help="The number of MCMC chains to use",
        type=int,
        default=default_chains,
    )

    parser.add_argument(
        "-i",
        "--iterations",
        help="""The number of MCMC iterations to use for
        each chain.""",
        type=int,
        default=default_iterations,
    )

    parser.add_argument(
        "-p",
        "--num-cpus",
        help="""The number of CPUs to use. Each read
        length will be processed in its own thread (so this is the maximum number of CPUs
        that is useful).""",
        type=int,
        default=1,
    )

    parser.add_argument("--type-field", default=default_type_field)

    parser.add_argument("--count-field", default=default_count_field)

    parser.add_argument("--position-field", default=default_position_field)

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    compile_rpbp_models.compile()

    # we will parallelize based on the lengths. So we need to know which lengths
    # are present in the metagene profiles file
    metagene_profiles = pd.read_csv(args.metagene_profiles)
    lengths = list(metagene_profiles["length"].unique())

    length_str = ",".join(str(int(l)) for l in lengths)
    msg = "Estimating Bayes factors for lengths: {}".format(length_str)
    logger.info(msg)

    length_groups = metagene_profiles.groupby("length")

    all_profile_estimates_df = parallel.apply_parallel_groups(
        length_groups,
        args.num_cpus,
        estimate_profile_bayes_factors,
        args,
        progress_bar=True,
    )

    msg = "Combining estimates into one data frame"
    logger.info(msg)

    all_profile_estimates_df = utils.remove_nones(all_profile_estimates_df)
    all_profile_estimates_df = pd.concat(all_profile_estimates_df)

    pandas_utils.write_df(all_profile_estimates_df, args.out, index=False)


if __name__ == "__main__":
    main()
