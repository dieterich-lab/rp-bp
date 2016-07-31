#! /usr/bin/env python3

import argparse
import logging
import os
import sys

import yaml
import misc.bio as bio
import misc.slurm as slurm
import misc.utils as utils

import riboutils.ribo_filenames as filenames

default_star_executable = "STAR"

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script creates all of the files necessary for downstream "
        "analysis performed with the rpbp package.")
    parser.add_argument('config', help="The (yaml) config file")

    parser.add_argument('--star-executable', help="The name of the STAR executable",
        default=default_star_executable)
    
    parser.add_argument('--overwrite', help="If this flag is present, existing files "
        "will be overwritten.", action='store_true')
    
    slurm.add_sbatch_options(parser)
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
                 'split-bed12-blocks',
                 args.star_executable
                ]
    utils.check_programs_exist(programs)

    
    required_keys = [   'genome_base_path',
                        'genome_name',
                        'gtf',
                        'fasta',
                        'ribosomal_fasta',
                        'ribosomal_index',
                        'star_index'
                    ]
    utils.check_keys_exist(config, required_keys)

    # now, check if we want to use slurm
    if args.use_slurm:
        cmd = ' '.join(sys.argv)
        slurm.check_sbatch(cmd, args=args)
        return
   
    # the rrna index
    cmd = "bowtie2-build-s {} {}".format(config['ribosomal_fasta'], config['ribosomal_index'])

    in_files = [config['ribosomal_fasta']]
    out_files = bio.get_bowtie2_index_files(config['ribosomal_index'])
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)
    
    # the STAR index
    mem = utils.human2bytes(args.mem)
    cmd = ("{} --runMode genomeGenerate --genomeDir {} --genomeFastaFiles {} "
        "--runThreadN {} --limitGenomeGenerateRAM {}".format(args.star_executable, 
        config['star_index'], config['fasta'], args.num_cpus, mem))
        
    in_files = [config['fasta']]
    out_files = bio.get_star_index_files(config['star_index'])
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)

    # extract the transcript fasta
    transcript_fasta = filenames.get_transcript_fasta(config['genome_base_path'], config['genome_name'])
    
    cmd = "gffread -W -w {} -g {} {}".format(transcript_fasta, config['fasta'], config['gtf'])
    in_files = [config['gtf'], config['fasta']]
    out_files = [transcript_fasta]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)
         
    # the deduplicated, annotated, genomic coordinates orfs
    orfs_genomic = filenames.get_orfs(config['genome_base_path'], config['genome_name'], 
        note=config.get('orf_note'))
    start_codons_str = utils.get_config_argument(config, 'start_codons')
    stop_codons_str = utils.get_config_argument(config, 'stop_codons')
    novel_id_str = utils.get_config_argument(config, 'novel_id_re')

    ignore_parsing_errors_str = ""
    if 'ignore_parsing_errors' in config:
        ignore_parsing_errors_str = "--ignore-parsing-errors"

    cmd = ("extract-orfs {} {} {} --num-cpus {} {} {} {} {}".format(transcript_fasta, orfs_genomic, 
        logging_str, args.num_cpus, start_codons_str, stop_codons_str, novel_id_str, ignore_parsing_errors_str))
    in_files = [transcript_fasta]
    out_files = [orfs_genomic]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)

    exons_file = filenames.get_exons(config['genome_base_path'], config['genome_name'],
        note=config.get('orf_note'))

    cmd = ("split-bed12-blocks {} {} --num-cpus {} {}".format(orfs_genomic, 
        exons_file, args.num_cpus, logging_str))
    in_files = [orfs_genomic]
    out_files = [exons_file]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)

if __name__ == '__main__':
    main()
