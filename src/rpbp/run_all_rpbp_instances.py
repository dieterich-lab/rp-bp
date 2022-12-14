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

import pbiotools.ribo.ribo_utils as ribo_utils

from rpbp.defaults import default_num_cpus, default_mem, star_executable

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""This is a helper script to submit a set of
        samples to SLURM. It can also be used to run a set of samples sequentially. Due to limitations
        on the config file specification, all of the samples must use the same reference indices
        obtained by running 'create-base-genome-profile.""",
    )

    parser.add_argument("config", help="The (yaml) configuration file")

    parser.add_argument("--tmp", help="The temp directory", default=None)

    parser.add_argument(
        "--overwrite",
        help="""If this flag is present, existing files
        will be overwritten.""",
        action="store_true",
    )

    parser.add_argument(
        "--profiles-only",
        help="""If this flag is present, then only
        the pre-processing part of the pipeline will be called, i.e. profiles
        will be created for each sample specified in the config file, but no predictions
        will be made.""",
        action="store_true",
    )

    parser.add_argument(
        "--merge-replicates",
        help="""If this flag is present, then
        the ORF profiles from the replicates will be merged before making the final
        predictions""",
        action="store_true",
    )

    parser.add_argument(
        "--run-replicates",
        help="""If this flag is given with the
        --merge-replicates flag, then both the replicates and the individual
        samples will be run. This flag has no effect if --merge-replicates is not
        given.""",
        action="store_true",
    )

    parser.add_argument(
        "-k",
        "--keep-intermediate-files",
        help="""If this flag is given,
        then all intermediate files will be kept; otherwise, they will be
        deleted. This feature is implemented piecemeal. If the --do-not-call flag
        is given, then nothing will be deleted. NOTE: intermediate files are
        necessary to perform read filtering QC!""",
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
