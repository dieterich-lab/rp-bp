#! /usr/bin/env python3

"""This is the main Rp-Bp script.

Calls:
    create-orf-profiles
    predict-translated-orfs (if running the prediction pipeline)
"""

import argparse
import logging
import shlex
import sys

import yaml

import pbiotools.utils.pgrm_utils as pgrm_utils
import pbiotools.misc.logging_utils as logging_utils
import pbiotools.misc.shell_utils as shell_utils
import pbiotools.misc.slurm as slurm
import pbiotools.misc.utils as utils

from rpbp.defaults import default_num_cpus, default_mem, star_executable

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""This script runs the Rp-Bp pipelines
        on a given sample. It requires a YAML config file that includes a number of keys.
        Please see the documentation for a complete description.""",
    )

    parser.add_argument("raw_data", help="The raw data file (fastq[.gz])")

    parser.add_argument("config", help="The (yaml) configuration file")

    parser.add_argument(
        "name", help="The name for the dataset, used in the created files"
    )

    parser.add_argument("--tmp", help="The temp directory", default=None)

    parser.add_argument(
        "--overwrite",
        help="If this flag is present, existing files " "will be overwritten.",
        action="store_true",
    )

    parser.add_argument(
        "--profiles-only",
        help="""If this flag is present, then only
        the ORF profiles will be created""",
        action="store_true",
    )

    parser.add_argument(
        "-k",
        "--keep-intermediate-files",
        help="""If this flag is given,
        then all intermediate files will be kept; otherwise, they will be
        deleted. This feature is implemented piecemeal. If the --do-not-call flag
        is given, then nothing will be deleted.""",
        action="store_true",
    )

    parser.add_argument(
        "--write-unfiltered",
        help="""If this flag is given, in addition to the default
        filtered predictions (longest ORF for each stop codon, then
        highest Bayes factor among overlapping ORFs), output all ORF
        predictions""",
        action="store_true",
    )

    slurm.add_sbatch_options(parser, num_cpus=default_num_cpus, mem=default_mem)
    logging_utils.add_logging_options(parser)
    pgrm_utils.add_star_options(parser, star_executable)
    pgrm_utils.add_flexbar_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    config = yaml.load(open(args.config), Loader=yaml.FullLoader)

    # check that all of the necessary programs are callable
    programs = [
        "flexbar",
        args.star_executable,
        "samtools",
        "bowtie2",
        "create-base-genome-profile",
        "remove-multimapping-reads",
        "extract-metagene-profiles",
        "estimate-metagene-profile-bayes-factors",
        "select-periodic-offsets",
        "extract-orf-profiles",
        "estimate-orf-bayes-factors",
        "select-final-prediction-set",
        "create-orf-profiles",
        "predict-translated-orfs",
    ]
    shell_utils.check_programs_exist(programs)

    required_keys = [
        "riboseq_data",
        "ribosomal_index",
        "star_index",
        "genome_base_path",
        "genome_name",
        "fasta",
        "gtf",
    ]
    utils.check_keys_exist(config, required_keys)

    # if using slurm, submit the script, but we cannot use sys.argv directly
    # as the shell strips the quotes around the arguments
    if args.use_slurm:
        cmd = "{}".format(" ".join("'" + s + "'" if '"' in s else s for s in sys.argv))
        slurm.check_sbatch(cmd, args=args)
        return

    # handle all option strings to call programs
    logging_str = logging_utils.get_logging_options_string(args)
    star_str = pgrm_utils.get_star_options_string(args)
    flexbar_str = pgrm_utils.get_flexbar_options_string(args)

    # handle do_not_call so that we do call the preprocessing script,
    # but that it does not run anything
    call = not args.do_not_call
    do_not_call_str = ""
    if not call:
        do_not_call_str = "--do-not-call"

    overwrite_str = ""
    if args.overwrite:
        overwrite_str = "--overwrite"

    keep_intermediate_str = ""
    if args.keep_intermediate_files:
        keep_intermediate_str = "--keep-intermediate-files"

    unfiltered_str = ""
    if args.write_unfiltered:
        unfiltered_str = "--write-unfiltered"

    tmp_str = ""
    if args.tmp is not None:
        tmp_str = "--tmp {}".format(shlex.quote(args.tmp))

    mem_str = "--mem {}".format(shlex.quote(args.mem))

    cmd = "create-orf-profiles {} {} {} --num-cpus {} {} {} {} {} {} {} {} {}".format(
        args.raw_data,
        args.config,
        args.name,
        args.num_cpus,
        mem_str,
        do_not_call_str,
        overwrite_str,
        keep_intermediate_str,
        logging_str,
        tmp_str,
        star_str,
        flexbar_str,
    )

    shell_utils.check_call(cmd)

    # check if we only want to create the profiles
    if args.profiles_only:
        return

    # then we predict the ORFs
    cmd = "predict-translated-orfs {} {} --num-cpus {} {} {} {} {}".format(
        args.config,
        args.name,
        args.num_cpus,
        do_not_call_str,
        overwrite_str,
        unfiltered_str,
        logging_str,
    )
    shell_utils.check_call(cmd)


if __name__ == "__main__":
    main()
