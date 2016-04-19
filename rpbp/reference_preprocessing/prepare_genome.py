#! /usr/bin/env python3

import argparse
import yaml
import os
import logging

import misc.bio as bio
import misc.utils as utils

import rpbp.filenames as filenames

default_num_procs = 2
default_mem = "8G"

default_star_executable = "STAR"
default_sjdb_overhang = 50

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
                 'extract-orfs',
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
    
    cmd = "gffread -W -w {} -g {} {}".format(transcript_fasta, config['fasta'], config['gtf'])
    in_files = [config['gtf'], config['fasta']]
    out_files = [transcript_fasta]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)
         
    # the deduplicated, annotated, genomic coordinates orfs
    orfs_genomic = filenames.get_orfs(config['base_path'], config['name'], note=config.get('orf_note'))
    start_codons_str = utils.get_config_argument(config, 'start_codons')
    stop_codons_str = utils.get_config_argument(config, 'stop_codons')
    novel_id_str = utils.get_config_argument(config, 'novel_id_re')

    ignore_parsing_errors_str = ""
    if 'ignore_parsing_errors' in config:
        ignore_parsing_errors_str = "--ignore-parsing-errors"

    cmd = ("extract-orfs {} {} {} --num-procs {} {} {} {} {}".format(transcript_fasta, orfs_genomic, 
        logging_str, args.num_procs, start_codons_str, stop_codons_str, novel_id_str, ignore_parsing_errors_str))
    in_files = [transcript_fasta]
    out_files = [orfs_genomic]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)

    # the rrna index
    in_files = [config['ribosomal_fasta']]
    out_files = bio.get_bowtie2_index_files(config['ribosomal_index'])
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)

    # the STAR index
    sjdb_overhang_str = utils.get_config_argument(config, 'sjdb_overhang', 'sjdbOverhang', 
        default=default_sjdb_overhang)
    mem = utils.human2bytes(args.mem)
    cmd = ("{} --runMode genomeGenerate --genomeDir {} --genomeFastaFiles {} --sjdbGTFfile {} "
        "{} --runThreadN {} --limitGenomeGenerateRAM {}".format(args.star_executable, 
        config['star_index'], config['fasta'], config['gtf'], sjdb_overhang_str, args.num_procs, mem))
        
    in_files = [config['gtf'], config['fasta']]
    out_files = bio.get_star_index_files(config['star_index'])
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)

if __name__ == '__main__':
    main()
