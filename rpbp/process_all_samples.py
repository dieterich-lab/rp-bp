#! /usr/bin/env python3

import argparse
import yaml
import misc.utils as utils

default_num_procs = 2
default_tmp = utils.abspath('tmp')
default_star_executable = "STAR"

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This is a helper script which submits a set of samples to SLURM. It "
        "can also be used to run a set of samples sequentially. Due to limitations on "
        "the config file specification, all of the samples must use the same reference "
        "indices (i.e., genome sequence, set of ORFs, etc.).")

    parser.add_argument('config', help="The (yaml) config file")

    parser.add_argument('-p', '--num-procs', help="The number of processors to use",
        type=int, default=default_num_procs)
    parser.add_argument('--tmp', help="The temp directory for pybedtools", default=default_tmp)
    parser.add_argument('--star-executable', help="The name of the STAR executable",
        default=default_star_executable)

    parser.add_argument('--do-not-call', action='store_true')
    parser.add_argument('--overwrite', help="If this flag is present, existing files "
        "will be overwritten.", action='store_true')
    
    utils.add_sbatch_options(parser)
    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    logging_str = utils.get_logging_options_string(args)

    config = yaml.load(open(args.config))
    call = not args.do_not_call


    # check that all of the necessary programs are callable
    programs =  [
                    'flexbar',
                    args.star_executable,
                    'samtools',
                    'bowtie2',
                    'bamToBed',
                    'fastaFromBed',
                    'create-base-genome-profile',
                    'remove-multimapping-reads',
                    'filter-alignment-lengths',
                    'extract-metagene-profiles',
                    'estimate-metagene-profile-bayes-factors',
                    'select-periodic-offsets',
                    'extract-orf-profiles',
                    'estimate-orf-bayes-factors',
                    'select-final-prediction-set',
                    'create-filtered-genome-profile',
                    'predict-translated-orfs',
                    'run-rpbp-pipeline'
                ]
    utils.check_programs_exist(programs)

    
    required_keys = [   
                        'riboseq_data',
                        'ribosomal_index',
                        'star_index',
                        'fasta',
                        'gtf',
                        'orfs',
                        'translated_models',
                        'untranslated_models',
                        'periodic_models',
                        'nonperiodic_models'
                    ]
    utils.check_keys_exist(config, required_keys)

    note_str = config.get('note', None)

    # the first step is the standard riboseq preprocessing
    
    # handle do_not_call so that we _do_ call the preprocessing script, but that it does not run anything
    do_not_call_str = ""
    if not call:
        do_not_call_str = "--do-not-call"

    overwrite_str = ""
    if args.overwrite:
        overwrite_str = "--overwrite"
    
    star_str = "--star-executable {}".format(args.star_executable)

    for name, data in config['samples']:

        cmd = "run-rpbp-pipeline {} {} {} --num-procs {} --tmp {} {} {} {} {}".format(data, 
                args.config, name, args.num_procs, args.tmp, do_not_call_str, 
                overwrite_str, logging_str, star_str)

        utils.check_sbatch(cmd, args=args)

if __name__ == '__main__':
    main()
