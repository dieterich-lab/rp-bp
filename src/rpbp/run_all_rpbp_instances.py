#! /usr/bin/env python3

import argparse
import logging
import os
import shlex
import yaml

from collections import defaultdict

import pbiotools.utils.pgrm_utils as pgrm_utils
import pbiotools.misc.logging_utils as logging_utils
import pbiotools.misc.shell_utils as shell_utils
import pbiotools.misc.slurm as slurm
import pbiotools.misc.utils as utils

import rpbp.ribo_utils.utils as ribo_utils

from rpbp.defaults import default_num_cpus, default_mem, star_executable

logger = logging.getLogger(__name__)


def get_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Call Rp-Bp for each sample in the configuration file.""",
    )

    parser.add_argument("config", help="A YAML configuration file.")

    parser.add_argument(
        "--profiles-only",
        help="Run the periodicity estimation only (ORF profile construction) "
        "for each sample in the configuration file.",
        action="store_true",
    )

    parser.add_argument(
        "--merge-replicates",
        help="Predict Ribo-seq ORFs in merged profiles. If ``--merge-replicates``, "
        "then use ``--run-replicates`` to also predict Ribo-seq ORFs in all samples.",
        action="store_true",
    )

    parser.add_argument(
        "--run-replicates",
        help="With ``--merge-replicates``, predict Ribo-seq ORFs in all samples and "
        "in merged profiles. This has no effect without ``--merge-replicates``, *i.e.* "
        "predictions are made for all samples by default.",
        action="store_true",
    )

    parser.add_argument(
        "--write-unfiltered",
        help="In addition to the default filtered predictions, output "
        "all overlapping predicted Ribo-seq ORFs.",
        action="store_true",
    )

    parser.add_argument(
        "-k",
        "--keep-intermediate-files",
        help="Keep all intermediate files. Intermediate files are "
        "necessary to perform read filtering quality control.",
        action="store_true",
    )

    parser.add_argument(
        "--overwrite",
        help="Overwrite existing output.",
        action="store_true",
    )

    parser.add_argument("--tmp", help="A temporary directory (STAR).", default=None)

    pgrm_utils.add_flexbar_options(parser)
    pgrm_utils.add_star_options(parser, star_executable)
    slurm.add_sbatch_options(parser, num_cpus=default_num_cpus, mem=default_mem)
    logging_utils.add_logging_options(parser)

    return parser


def main():
    parser = get_parser()
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
        "run-rpbp-pipeline",
    ]
    shell_utils.check_programs_exist(programs)

    required_keys = [
        "riboseq_data",
        "riboseq_samples",
        "ribosomal_index",
        "star_index",
        "genome_base_path",
        "genome_name",
        "fasta",
        "gtf",
    ]
    utils.check_keys_exist(config, required_keys)

    # handle all option strings to call the pipeline script
    logging_str = logging_utils.get_logging_options_string(args)
    star_str = pgrm_utils.get_star_options_string(args)
    flexbar_str = pgrm_utils.get_flexbar_options_string(args)

    # handle do_not_call so that we do call the pipeline script, but that it does not run anything
    call = not args.do_not_call
    do_not_call_str = ""
    if not call:
        do_not_call_str = "--do-not-call"
    args.do_not_call = False

    overwrite_str = ""
    if args.overwrite:
        overwrite_str = "--overwrite"

    mem_str = "--mem {}".format(shlex.quote(args.mem))

    keep_intermediate_str = ""
    if args.keep_intermediate_files:
        keep_intermediate_str = "--keep-intermediate-files"

    unfiltered_str = ""
    if args.write_unfiltered:
        unfiltered_str = "--write-unfiltered"

    # check if we only want to create the profiles, in this case
    # we call run-rpbp-pipeline with the --profiles-only option
    profiles_only_str = ""
    if args.profiles_only:
        if args.merge_replicates:
            msg = (
                "The --profiles-only option was given, this option has"
                "precedence, and it will override the --merge-replicates option!"
            )
            logger.warning(msg)
        args.merge_replicates = False
        profiles_only_str = "--profiles-only"

    # if we merge the replicates, then we only use the rpbp script to create
    # the ORF profiles, but we still make predictions
    if args.merge_replicates and not args.run_replicates:
        profiles_only_str = "--profiles-only"

    if args.run_replicates and not args.merge_replicates:
        msg = (
            "The --run-replicates option was given without the --merge-replicates "
            "option. It will be ignored."
        )
        logger.warning(msg)

    # collect the job_ids in case we are using slurm and need to merge replicates
    rep_to_condition = ribo_utils.get_riboseq_replicates_reverse_map(config)
    job_ids_mapping = defaultdict(list)

    sample_names = sorted(config["riboseq_samples"].keys())

    for sample_name in sample_names:
        data = config["riboseq_samples"][sample_name]

        tmp_str = ""
        if args.tmp is not None:
            tmp = os.path.join(args.tmp, "{}_rpbp".format(sample_name))
            tmp_str = "--tmp {}".format(tmp)

        cmd = "run-rpbp-pipeline {} {} {} --num-cpus {} {} {} {} {} {} {} {} {} {} {}".format(
            data,
            args.config,
            sample_name,
            args.num_cpus,
            mem_str,
            tmp_str,
            do_not_call_str,
            overwrite_str,
            profiles_only_str,
            keep_intermediate_str,
            unfiltered_str,
            logging_str,
            star_str,
            flexbar_str,
        )

        job_id = slurm.check_sbatch(cmd, args=args)
        job_ids_mapping[rep_to_condition[sample_name]].append(job_id)

    # now, if we are running the "standard" pipeline, we are done
    if not args.merge_replicates:
        return

    # otherwise, we need to merge the replicates for each condition
    riboseq_replicates = ribo_utils.get_riboseq_replicates(config)
    merge_replicates_str = "--merge-replicates"

    for condition_name in sorted(riboseq_replicates.keys()):

        # then we predict the ORFs
        cmd = "predict-translated-orfs {} {} --num-cpus {} {} {} {} {} {}".format(
            args.config,
            condition_name,
            args.num_cpus,
            do_not_call_str,
            overwrite_str,
            unfiltered_str,
            logging_str,
            merge_replicates_str,
        )

        job_ids = job_ids_mapping[condition_name]
        slurm.check_sbatch(cmd, args=args, dependencies=job_ids)


if __name__ == "__main__":
    main()
