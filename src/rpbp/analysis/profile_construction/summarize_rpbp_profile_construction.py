#! /usr/bin/env python3

import os
import sys
import argparse
import logging
import yaml
import csv

from pathlib import Path

import scipy.stats

import numpy as np
import pandas as pd

import pbiotools.misc.logging_utils as logging_utils
import pbiotools.misc.pandas_utils as pandas_utils
import pbiotools.misc.parallel as parallel
import pbiotools.misc.shell_utils as shell_utils
import pbiotools.misc.utils as utils

import rpbp.ribo_utils.utils as ribo_utils
import rpbp.ribo_utils.filenames as filenames

from rpbp.defaults import default_num_cpus, metagene_options

logger = logging.getLogger(__name__)


def get_read_filtering_summary(filename, args):

    overwrite_str = ""
    if args.overwrite:
        overwrite_str = "--overwrite"

    logging_str = logging_utils.get_logging_options_string(args)

    cpus_str = "--num-cpus {}".format(args.num_cpus)
    cmd = "get-all-read-filtering-counts {} {} {} {} {}".format(
        args.config, filename, overwrite_str, cpus_str, logging_str
    )
    in_files = [args.config]
    out_files = [filename]
    shell_utils.call_if_not_exists(
        cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=True
    )


def get_read_length_distributions(sample, config, is_unique, note, args):

    logging_str = logging_utils.get_logging_options_string(args)
    cpus_str = "--num-cpus {}".format(args.num_cpus)

    # distribution for this sample
    read_length_distribution = filenames.get_riboseq_read_length_distribution(
        config["riboseq_data"], sample, note=note
    )

    #  all aligned reads
    genome_bam = filenames.get_riboseq_bam(
        config["riboseq_data"], sample, is_unique=False, note=note
    )
    if is_unique:
        # uniquely aligned reads
        unique_bam = filenames.get_riboseq_bam(
            config["riboseq_data"], sample, is_unique=is_unique, note=note
        )
        in_files = [genome_bam, unique_bam]
        cmd = "get-read-length-distribution {} {} --out {} {} {}".format(
            genome_bam, unique_bam, read_length_distribution, logging_str, cpus_str
        )
    else:
        in_files = [genome_bam]
        cmd = "get-read-length-distribution {} --out {} {} {}".format(
            genome_bam, read_length_distribution, logging_str, cpus_str
        )

    out_files = [read_length_distribution]
    shell_utils.call_if_not_exists(
        cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=True
    )


def collect_all_read_length_distributions(sample, config, note, args):

    # distribution for this sample
    read_length_distribution = filenames.get_riboseq_read_length_distribution(
        config["riboseq_data"], sample, note=note
    )
    read_length_distribution = pd.read_csv(read_length_distribution)
    ret = pd.pivot_table(
        read_length_distribution, values="count", columns="length", index="basename"
    )
    return ret


def collect_lengths_and_offsets(sample, config, is_unique, note, args):

    # not only the periodic lengths and offsets, but all
    # and add status information based on filters

    min_metagene_profile_count = config.get(
        "min_metagene_profile_count", metagene_options["min_metagene_profile_count"]
    )
    min_bf_mean = config.get(
        "min_metagene_bf_mean", metagene_options["min_metagene_bf_mean"]
    )
    max_bf_var = config.get(
        "max_metagene_bf_var", metagene_options["max_metagene_bf_var"]
    )
    min_bf_likelihood = config.get(
        "min_metagene_bf_likelihood", metagene_options["min_metagene_bf_likelihood"]
    )

    periodic_offsets = filenames.get_periodic_offsets(
        config["riboseq_data"], sample, is_unique=is_unique, note=note
    )
    offsets_df = pd.read_csv(periodic_offsets)

    min_read_length = int(offsets_df["length"].min())
    max_read_length = int(offsets_df["length"].max())

    lengths, offsets, statuses = [], [], []
    for length in range(min_read_length, max_read_length + 1):

        # check which offset is used
        # select the row for this length
        mask_length = offsets_df["length"] == length

        # TODO: this is sometimes length 0. why?
        if sum(mask_length) == 0:
            continue

        length_row = offsets_df[mask_length].iloc[0]

        # now, check all of the filters
        offset = int(length_row["highest_peak_offset"])
        offset_status = "Used for analysis"

        if max_bf_var is not None:
            if (length_row["highest_peak_bf_mean"] <= min_bf_mean) or length_row[
                "highest_peak_bf_var"
            ] >= max_bf_var:
                offset_status = "BF mean too small or BF var too high"

        if min_bf_likelihood is not None:
            likelihood = 1 - scipy.stats.norm.cdf(
                min_bf_mean,
                length_row["highest_peak_bf_mean"],
                np.sqrt(length_row["highest_peak_bf_var"]),
            )
            if likelihood <= min_bf_likelihood:
                offset_status = "Likehood too small"

        if (max_bf_var is None) and (min_bf_likelihood is None):
            if length_row["highest_peak_bf_mean"] <= min_bf_mean:
                offset_status = "BF mean too small"

        if length_row["highest_peak_profile_sum"] < min_metagene_profile_count:
            offset_status = "Count too small"

        lengths.append(length)
        offsets.append(offset)
        statuses.append(offset_status)

    data = {"sample": sample, "length": lengths, "offset": offsets, "status": statuses}

    return pd.DataFrame(data, columns=["sample", "length", "offset", "status"])


def get_profile(sample, config, is_unique, note):
    """Get the name of the profile file from the given parameters."""

    lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(
        config, sample, is_unique=is_unique, default_params=metagene_options
    )

    if len(lengths) == 0:
        msg = "No periodic read lengths and offsets were found!"
        logger.critical(msg)
        return

    profiles = filenames.get_riboseq_profiles(
        config["riboseq_data"],
        sample,
        length=lengths,
        offset=offsets,
        is_unique=is_unique,
        note=note,
    )

    return profiles


def get_frame_counts(sample, config, is_unique, note):

    msg = "{}: extracting frame counts".format(sample)
    logger.info(msg)

    mtx = get_profile(sample, config, is_unique, note)

    # we don't need to load as sparse matrix, use numpy and
    # mask based on ORF offset (2nd column), taking into account that
    # mtx format is 1-based
    # skip header, including mtx format specifications

    profiles = np.loadtxt(mtx, skiprows=3)

    m_frame1 = np.where(profiles[:, 1] % 3 == 1)
    m_frame2 = np.where(profiles[:, 1] % 3 == 2)
    m_frame3 = np.where(profiles[:, 1] % 3 == 0)

    frame1 = profiles[m_frame1][:, 2].sum()
    frame2 = profiles[m_frame2][:, 2].sum()
    frame3 = profiles[m_frame3][:, 2].sum()

    ret = {"sample": sample, "frame": frame1, "frame+1": frame2, "frame+2": frame3}

    return pd.Series(ret)


def create_fastqc_reports(sample_data, config, note, is_unique, args):

    sample, raw_data = sample_data
    msg = "{}: creating fastqc reports".format(sample)
    logger.info(msg)

    # first get the filenames
    without_adapters = filenames.get_without_adapters_fastq(
        config["riboseq_data"], sample, note=note
    )
    with_rrna = filenames.get_with_rrna_fastq(config["riboseq_data"], sample, note=note)
    without_rrna = filenames.get_without_rrna_fastq(
        config["riboseq_data"], sample, note=note
    )
    genome_bam = filenames.get_riboseq_bam(config["riboseq_data"], sample, note=note)
    unique_bam = filenames.get_riboseq_bam(
        config["riboseq_data"], sample, is_unique=is_unique, note=note
    )

    raw_data_fastqc = filenames.get_raw_data_fastqc_data(
        config["riboseq_data"], raw_data
    )
    without_adapters_fastqc = filenames.get_without_adapters_fastqc_data(
        config["riboseq_data"], sample, note=note
    )
    with_rrna_fastqc = filenames.get_with_rrna_fastqc_data(
        config["riboseq_data"], sample, note=note
    )
    without_rrna_fastqc = filenames.get_without_rrna_fastqc_data(
        config["riboseq_data"], sample, note=note
    )

    genome_bam_fastqc = filenames.get_riboseq_bam_fastqc_data(
        config["riboseq_data"], sample, note=note
    )
    unique_bam_fastqc = filenames.get_riboseq_bam_fastqc_data(
        config["riboseq_data"], sample, is_unique=is_unique, note=note
    )

    # create the fastqc reports if they do not exist
    raw_data_fastqc_path = filenames.get_raw_data_fastqc_path(config["riboseq_data"])
    without_adapters_fastqc_path = filenames.get_without_adapters_fastqc(
        config["riboseq_data"]
    )
    with_rrna_fastqc_path = filenames.get_with_rrna_fastqc(config["riboseq_data"])
    without_rrna_fastqc_path = filenames.get_without_rrna_fastqc(config["riboseq_data"])
    without_rrna_mapping_fastqc_path = filenames.get_riboseq_bam_fastqc_path(
        config["riboseq_data"]
    )

    fastqc_tmp_str = ""
    if args.tmp is not None:
        fastqc_tmp_str = "--dir {}".format(args.tmp)

    all_infiles = [
        raw_data,
        without_adapters,
        with_rrna,
        without_rrna,
        genome_bam,
        unique_bam,
    ]
    all_fastqc_reports = [
        raw_data_fastqc,
        without_adapters_fastqc,
        with_rrna_fastqc,
        without_rrna_fastqc,
        genome_bam_fastqc,
        unique_bam_fastqc,
    ]
    all_paths = [
        raw_data_fastqc_path,
        without_adapters_fastqc_path,
        with_rrna_fastqc_path,
        without_rrna_fastqc_path,
        without_rrna_mapping_fastqc_path,
        without_rrna_mapping_fastqc_path,
    ]

    all_files = zip(all_infiles, all_fastqc_reports, all_paths)
    for inf, outf, fastqc_data in all_files:
        cmd = "fastqc --outdir {} --extract {} {}".format(
            fastqc_data, inf, fastqc_tmp_str
        )
        in_files = [inf]
        out_files = [outf]
        shell_utils.call_if_not_exists(
            cmd, out_files, in_files=in_files, overwrite=args.overwrite
        )

    # in some cases, fastqc can fail: make sure all of the reports are present
    missing_files = [f for f in all_fastqc_reports if not os.path.exists(f)]
    if len(missing_files) > 0:
        msg = "The following fastqc reports were not created correctly:\n"
        msg += "\n".join(missing_files)
        logger.warning(msg)


def get_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Summarize the ORF profile creation step "
        "and prepare data for the web application.",
    )

    parser.add_argument(
        "config",
        help="A YAML configuration file. " "The same used to run the pipeline.",
    )

    # TODO MultiQC
    parser.add_argument(
        "-c",
        "--create-fastqc-reports",
        help="Create FastQC reports.",
        action="store_true",
    )

    parser.add_argument(
        "--overwrite",
        help="Overwrite existing output",
        action="store_true",
    )

    parser.add_argument(
        "--tmp",
        help="A temporary directory (FastQC).",
        default=None,
    )

    parser.add_argument(
        "-p",
        "--num-cpus",
        help="The number of CPUs to use.",
        type=int,
        default=default_num_cpus,
    )

    logging_utils.add_logging_options(parser)

    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "[summarize-rpbp-profile-construction]: {}".format(" ".join(sys.argv))
    logger.info(msg)

    config = yaml.load(open(args.config), Loader=yaml.FullLoader)

    # calls TODO check
    programs = [
        "get-all-read-filtering-counts",
        "samtools",
        "get-read-length-distribution",
    ]

    if args.create_fastqc_reports:
        programs.extend(["fastqc", "java"])

    shell_utils.check_programs_exist(programs)

    required_keys = ["riboseq_data", "riboseq_samples"]
    utils.check_keys_exist(config, required_keys)

    # create directory for summary data
    sub_folder = Path("analysis", "profile_construction")
    Path(config["riboseq_data"], sub_folder).mkdir(parents=True, exist_ok=True)

    # nomenclature
    project = config.get("project_name", "rpbp")
    note = config.get("note", None)
    is_unique = not config.get("keep_riboseq_multimappers", False)

    msg = "Collecting all read filtering counts..."
    logger.info(msg)

    # 1. read filtering
    read_filtering_counts = filenames.get_riboseq_read_filtering_counts(
        config["riboseq_data"], project, sub_folder=sub_folder.as_posix(), note=note
    )
    get_read_filtering_summary(read_filtering_counts, args)

    msg = "Collecting all read length distributions..."
    logger.info(msg)

    # 2. read length distributions
    samples = sorted(config["riboseq_samples"].keys())
    for sample in samples:
        get_read_length_distributions(sample, config, is_unique, note, args)
    # collect in one dataframe for the app
    all_length_distributions = parallel.apply_iter_simple(
        samples, collect_all_read_length_distributions, config, note, args
    )
    all_length_distributions = pd.concat(all_length_distributions).fillna(0)
    all_length_distributions["is_unique"] = "False"
    all_length_distributions.loc[
        all_length_distributions.index.str.contains("unique"), "is_unique"
    ] = "True"
    # adjust name - the function called gets the name from the file name
    all_length_distributions.index = all_length_distributions.index.str.replace(
        "-unique", ""
    )
    if note is not None:
        repl = ".{}".format(note)
        all_length_distributions.index = all_length_distributions.index.str.replace(
            repl, ""
        )
    all_length_distributions.index.name = "sample"
    # make prettier
    column = all_length_distributions.pop("is_unique")
    all_length_distributions.insert(0, "is_unique", column)

    summary_file = filenames.get_riboseq_read_length_distribution(
        config["riboseq_data"], project, sub_folder=sub_folder.as_posix(), note=note
    )
    pandas_utils.write_df(
        all_length_distributions,
        summary_file,
        index=True,
        sep=",",
        header=True,
        do_not_compress=False,
        quoting=csv.QUOTE_NONE,
    )

    msg = "Collecting all lengths and offsets..."
    logger.info(msg)

    # 3. collect all (periodic) lengths and offsets
    all_lengths_and_offsets = parallel.apply_iter_simple(
        samples, collect_lengths_and_offsets, config, is_unique, note, args
    )
    all_lengths_and_offsets = pd.concat(all_lengths_and_offsets)

    summary_file = filenames.get_periodic_offsets(
        config["riboseq_data"],
        project,
        sub_folder=sub_folder.as_posix(),
        is_unique=is_unique,
        note=note,
    )
    pandas_utils.write_df(
        all_lengths_and_offsets,
        summary_file,
        index=False,
        sep=",",
        header=True,
        do_not_compress=False,
        quoting=csv.QUOTE_NONE,
    )

    msg = "Collecting all ORF profile counts per frame..."
    logger.info(msg)

    # 4. collect reads per frame using the ORF profiles - currently all ORFs
    frame_counts = parallel.apply_parallel_iter(
        samples,
        args.num_cpus,
        get_frame_counts,
        config,
        is_unique,
        note,
        progress_bar=True,
        backend="multiprocessing",
    )
    frame_counts = pd.DataFrame(frame_counts)

    summary_file = filenames.get_riboseq_frame_counts(
        config["riboseq_data"],
        project,
        sub_folder=sub_folder.as_posix(),
        is_unique=is_unique,
        note=note,
    )
    pandas_utils.write_df(
        frame_counts,
        summary_file,
        index=False,
        sep=",",
        header=True,
        do_not_compress=False,
        quoting=csv.QUOTE_NONE,
    )

    if args.create_fastqc_reports:

        msg = "Calling FastQC..."
        logger.info(msg)

        parallel.apply_parallel_iter(
            config["riboseq_samples"].items(),
            args.num_cpus,
            create_fastqc_reports,
            config,
            note,
            is_unique,
            args,
        )


if __name__ == "__main__":
    main()
