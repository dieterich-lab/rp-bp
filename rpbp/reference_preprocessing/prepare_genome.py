#! /usr/bin/env python3

import argparse
import yaml
import os
import logging

import misc.utils as utils

import rpbp.filenames as filenames

default_num_procs = 2
default_mem = "8G"

default_star_executable = "STAR"

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script creates all of the files necessary for downstream "
        "analysis performed with the rpbp package.")
    parser.add_argument('config', help="The (yaml) config file")

    parser.add_argument('-p', '--num-procs', help="The number of processors to use",
        type=int, default=default_num_procs)
    parser.add_argument('-m', '--mem', help="The amount of memory to use for creating "
        "the STAR index", default=default_mem)

    parser.add_argument('--star-executable', help="The name of the STAR executable",
        default=default_star_executable)
    
    parser.add_argument('--overwrite', help="If this flag is present, existing files "
        "will be overwritten.", action='store_true')
    parser.add_argument('--do-not-call', action='store_true')
    
    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    logging_str = utils.get_logging_options_string(args)

    config = yaml.load(open(args.config))
    call = not args.do_not_call

    # check that all of the necessary programs are callable
    programs =  [
                 'extract-transcript-fasta',
                 'extract-orfs',
                 'create-star-reference',
                 'bowtie2-build-s',
                 'gffread',
                 'intersectBed',
                 args.star_executable
                ]
    utils.check_programs_exist(programs)

    
    required_keys = [   'base_path',
                        'name',
                        'gtf',
                        'fasta',
                        'ribosomal_fasta',
                        'ribosomal_index',
                        'star_index'
                    ]
    utils.check_keys_exist(config, required_keys)


    # first, extract the transcript fasta
    transcript_fasta = filenames.get_transcript_fasta(config['base_path'], config['name'])
    cmd  = "extract-transcript-fasta {} {} {}".format(config['gtf'], config['fasta'], transcript_fasta)
    
    in_files = [config['gtf'], config['fasta']]
    out_files = [transcript_fasta]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)
         
    # the deduplicated, annotated, genomic coordinates orfs
    orfs_genomic = filenames.get_orfs(config['base_path'], config['name'], note=config.get('orf_note'))
    start_codons_str = utils.get_config_argument(config, 'start_codons')
    stop_codons_str = utils.get_config_argument(config, 'stop_codons')
    
    cmd = ("extract-orfs {} {} {} --num-procs {} {} {}".format(transcript_fasta, orfs_genomic, 
        logging_str, args.num_procs, start_codons_str, stop_codons_str))
    in_files = [transcript_fasta]
    out_files = [orfs_genomic]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)

    # the rrna index
    if not args.overwrite:
        msg = ("The --overwrite flag was not given, but the ribosomal bowtie2 "
            "index will still be re-created.")
        logging.warn(msg)

    cmd = "bowtie2-build-s {} {}".format(config['ribosomal_fasta'], config['ribosomal_index'])
    utils.check_call(cmd, call=call)

    # the STAR index
    if not args.overwrite:
        msg = ("The --overwrite flag was not given, but the STAR genome index "
            "will still be re-created.")
        logging.warn(msg)

    cmd  = "create-star-reference {} {} {} --num-procs {} --mem {} --star-executable {}".format(config['gtf'], 
        config['fasta'], config['star_index'], args.num_procs, args.mem, args.star_executable)
    utils.check_call(cmd, call=call)

if __name__ == '__main__':
    main()
