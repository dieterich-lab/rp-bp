#! /usr/bin/env python3

import argparse
import logging
import os
import shlex
import yaml

import pbio.utils.star_utils as star_utils
import pbio.misc.logging_utils as logging_utils
import pbio.misc.shell_utils as shell_utils
import pbio.misc.slurm as slurm
import pbio.misc.utils as utils

import pbio.ribo.ribo_utils as ribo_utils

logger = logging.getLogger(__name__)

default_num_procs = 1
default_tmp = None  # utils.abspath('tmp')


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="This is a helper script which submits a set of "
                                                 "samples to SLURM. It can also be used to run a "
                                                 "set of samples sequentially. Due to limitations on "
                                                 "the config file specification, all of the samples "
                                                 "must use the same reference indices (i.e., genome "
                                                 "sequence, set of ORFs, etc.).")

    parser.add_argument('config', help="The (yaml) config file")

    parser.add_argument('--tmp', help="The temp directory", default=default_tmp)

    parser.add_argument('--flexbar-options', help="""Optional argument: a space-delimited 
        list of options to pass to flexbar. Each option must be quoted separately as in 
        "--flexbarOption value", using hard, then soft quotes, where "--flexbarOption" 
        is the long parameter name from flexbar and "value" is the value given to this parameter. 
        If specified, flexbar options will override default settings.""", nargs='*', type=str)

    parser.add_argument('--overwrite', help="If this flag is present, existing files "
                                            "will be overwritten.", action='store_true')

    parser.add_argument('--profiles-only', help="""If this flag is present, then only
        the pre-processing part of the pipeline will be called, i.e. profiles
        will be created for each sample specified in the config file, but no predictions
        will be made.""", action='store_true')

    parser.add_argument('--merge-replicates', help="""If this flag is present, then
        the ORF profiles from the replicates will be merged before making the final
        predictions""", action='store_true')

    parser.add_argument('--run-replicates', help="""If this flag is given with the
        --merge-replicates flag, then both the replicates *and* the individual
        samples will be run. This flag has no effect if --merge-replicates is not
        given.""", action='store_true')
         
    parser.add_argument('-k', '--keep-intermediate-files', help="""If this flag is given,
        then all intermediate files will be kept; otherwise, they will be
        deleted. This feature is implemented piecemeal. If the --do-not-call flag
        is given, then nothing will be deleted.""", action='store_true')
    
    star_utils.add_star_options(parser)
    slurm.add_sbatch_options(parser)
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    logging_str = logging_utils.get_logging_options_string(args)
    star_str = star_utils.get_star_options_string(args)

    config = yaml.load(open(args.config), Loader=yaml.FullLoader)
    call = not args.do_not_call

    # check that all of the necessary programs are callable
    programs = [
        'flexbar',
        args.star_executable,
        'samtools',
        'bowtie2',
        'create-base-genome-profile',
        'remove-multimapping-reads',
        'extract-metagene-profiles',
        'estimate-metagene-profile-bayes-factors',
        'select-periodic-offsets',
        'extract-orf-profiles',
        'estimate-orf-bayes-factors',
        'select-final-prediction-set',
        'create-orf-profiles',
        'predict-translated-orfs',
        'run-rpbp-pipeline'
    ]
    shell_utils.check_programs_exist(programs)
    
    required_keys = [
        'riboseq_data',
        'riboseq_samples',
        'ribosomal_index',
        'star_index',
        'genome_base_path',
        'genome_name',
        'fasta',
        'gtf'
    ]
    utils.check_keys_exist(config, required_keys)

    note = config.get('note', None)

    # handle do_not_call so that we _do_ call the pipeline script, but that it does not run anything
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

    # check if we only want to create the profiles, in this case
    # we call run-rpbp-pipeline with the --profiles-only option
    profiles_only_str = ""
    if args.profiles_only:
        args.merge_replicates = False
        profiles_only_str = "--profiles-only"
        msg = ("The --profiles-only option was given, this will override --merge-replicates "
               "and/or --run-replicates, if these options were also given!")
        logger.info(msg)

    # if we merge the replicates, then we only use the rpbp script to create
    # the ORF profiles, but we still make predictions
    if args.merge_replicates and not args.run_replicates:
        profiles_only_str = "--profiles-only"

    if args.run_replicates and not args.merge_replicates:
        msg = ("The --run-replicates option was given without the --merge-replicates "
               "option. It will be ignored.")
        logger.warning(msg)
    
    tmp_str = ""
    if args.tmp is not None:
        tmp_str = "--tmp {}".format(args.tmp)

    flexbar_option_str = ""
    if args.flexbar_options is not None:
        flexbar_option_str = "--flexbar-options {}".format(' '.join("'" + flx_op + "'"
                                                                    for flx_op in args.flexbar_options))

    # collect the job_ids in case we are using slurm and need to merge replicates
    job_ids = []

    sample_names = sorted(config['riboseq_samples'].keys())

    for sample_name in sample_names:
        data = config['riboseq_samples'][sample_name]

        tmp_str = ""
        if args.tmp is not None:
            tmp = os.path.join(args.tmp, "{}_{}_rpbp".format(sample_name, note))
            tmp_str = "--tmp {}".format(tmp)

        cmd = "run-rpbp-pipeline {} {} {} --num-cpus {} {} {} {} {} {} {} {} {} {}".format(
            data, 
            args.config, 
            sample_name, 
            args.num_cpus, 
            tmp_str, 
            do_not_call_str, 
            overwrite_str, 
            logging_str, 
            star_str, 
            profiles_only_str,
            flexbar_option_str,
            keep_intermediate_str,
            mem_str
        )

        job_id = slurm.check_sbatch(cmd, args=args)

        job_ids.append(job_id)

    # now, if we are running the "standard" pipeline, we are finished
    if not args.merge_replicates:
        return

    # otherwise, we need to merge the replicates for each condition
    riboseq_replicates = ribo_utils.get_riboseq_replicates(config)
    merge_replicates_str = "--merge-replicates"

    for condition_name in sorted(riboseq_replicates.keys()):
    
        # then we predict the ORFs
        cmd = "predict-translated-orfs {} {} --num-cpus {} {} {} {} {}".format(
            args.config, 
            condition_name, 
            args.num_cpus, 
            do_not_call_str, 
            overwrite_str, 
            logging_str, 
            merge_replicates_str
        )

        slurm.check_sbatch(cmd, args=args, dependencies=job_ids)


if __name__ == '__main__':
    main()
