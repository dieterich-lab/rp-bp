#! /usr/bin/env python3

import argparse
import logging
import os
import sys

import yaml
import misc.bio as bio
import misc.bio_utils.bed_utils as bed_utils
import misc.logging_utils as logging_utils
import misc.shell_utils as shell_utils
import misc.slurm as slurm
import misc.utils as utils

import riboutils.ribo_filenames as filenames

logger = logging.getLogger(__name__)

default_star_executable = "STAR"

def get_orfs(gtf, args, config, is_annotated=False, is_de_novo=False):
    """ This helper function processes a GTF file into its ORFs.
    """
    call = not args.do_not_call
    chr_name_file = os.path.join(config['star_index'], 'chrName.txt')
    chr_name_str = "--chr-name-file {}".format(chr_name_file)

    logging_str = logging_utils.get_logging_options_string(args)
    cpus_str = "--num-cpus {}".format(args.num_cpus)

    # extract a bed12 of the annotated ORFs
    transcript_bed = filenames.get_bed(config['genome_base_path'], 
        config['genome_name'], is_merged=False, is_annotated=is_annotated, 
        is_de_novo=is_de_novo)
    
    cmd = ("gtf-to-bed12 {} {} --num-cpus {} {} {}".format(gtf, 
        transcript_bed, args.num_cpus, chr_name_str, logging_str))
    in_files = [gtf]
    out_files = [transcript_bed]
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files, 
        overwrite=args.overwrite, call=call)

    exons_file = filenames.get_exons(config['genome_base_path'], 
        config['genome_name'], is_annotated=is_annotated, is_de_novo=is_de_novo)

    cmd = ("split-bed12-blocks {} {} --num-cpus {} {}".format(transcript_bed, 
        exons_file, args.num_cpus, logging_str))
    in_files = [transcript_bed]
    out_files = [exons_file]
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files, 
        overwrite=args.overwrite, call=call)

    # extract the transcript fasta
    transcript_fasta = filenames.get_transcript_fasta(config['genome_base_path'], 
        config['genome_name'], is_annotated=is_annotated, is_de_novo=is_de_novo)

    cmd = ("extract-bed-sequences {} {} {} {}".format(transcript_bed, 
        config['fasta'], transcript_fasta, logging_str))
    in_files = [transcript_bed, config['fasta']]
    out_files = [transcript_fasta]
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files, 
        overwrite=args.overwrite, call=call)

    # new approach for extracting orfs
    orfs_genomic = filenames.get_orfs(config['genome_base_path'], config['genome_name'], 
        note=config.get('orf_note'), is_annotated=is_annotated, is_de_novo=is_de_novo)
    start_codons_str = utils.get_config_argument(config, 'start_codons')
    stop_codons_str = utils.get_config_argument(config, 'stop_codons')

    cmd = "extract-orf-coordinates {} {} {} {} {} {} {}".format(transcript_bed, 
        transcript_fasta, orfs_genomic, cpus_str, start_codons_str, 
        stop_codons_str, logging_str)

    in_files = [transcript_fasta, transcript_bed]
    out_files = [orfs_genomic]
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files, 
        overwrite=args.overwrite, call=call)

    exons_file = filenames.get_exons(config['genome_base_path'], config['genome_name'],
        note=config.get('orf_note'), is_annotated=is_annotated, is_de_novo=is_de_novo)

    cmd = ("split-bed12-blocks {} {} --num-cpus {} {}".format(orfs_genomic, 
        exons_file, args.num_cpus, logging_str))
    in_files = [orfs_genomic]
    out_files = [exons_file]
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files, 
        overwrite=args.overwrite, call=call)

    # label the orfs
    labeled_orfs = orfs_genomic # no need to keep the unannotated ones around

    # we always label wrt the annotated annotations
    annotated_bed = filenames.get_bed(config['genome_base_path'], 
        config['genome_name'], is_merged=False, is_annotated=True)

    de_novo_str = ""
    if is_de_novo:
         de_novo_str = "--label-prefix \"novel_\" --filter --nonoverlapping-label \"novel\""
    
    cmd = "label-orfs {} {} {} {} {} {} {}".format(annotated_bed, orfs_genomic, 
        exons_file, labeled_orfs, cpus_str, de_novo_str, logging_str)
    in_files = [annotated_bed, orfs_genomic, exons_file]
    #  since we are reusing the name, it will already exist
    out_files  = None # [] # [labeled_orfs]
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files, 
        overwrite=args.overwrite, call=call)



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
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    logging_str = logging_utils.get_logging_options_string(args)

    config = yaml.load(open(args.config))
    call = not args.do_not_call

    # check that all of the necessary programs are callable
    programs =  [
                 'extract-orfs',
                 'bowtie2-build-s',
                 'split-bed12-blocks',
                 'gtf-to-bed12',
                 args.star_executable
                ]
    shell_utils.check_programs_exist(programs)

    
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
    cmd = "bowtie2-build-s {} {}".format(config['ribosomal_fasta'], 
        config['ribosomal_index'])

    in_files = [config['ribosomal_fasta']]
    out_files = bio.get_bowtie2_index_files(config['ribosomal_index'])
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files, 
        overwrite=args.overwrite, call=call)
    
    # the STAR index
    mem = utils.human2bytes(args.mem)
    cmd = ("{} --runMode genomeGenerate --genomeDir {} --genomeFastaFiles {} "
        "--runThreadN {} --limitGenomeGenerateRAM {}".format(args.star_executable, 
        config['star_index'], config['fasta'], args.num_cpus, mem))
        
    in_files = [config['fasta']]
    out_files = bio.get_star_index_files(config['star_index'])
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files, 
        overwrite=args.overwrite, call=call)

    # get the main orfs
    get_orfs(config['gtf'], args, config, is_annotated=True, is_de_novo=False)

    # eventually, we will use these names
    annotated_orfs = filenames.get_orfs(config['genome_base_path'], 
        config['genome_name'], note=config.get('orf_note'), is_annotated=True,
        is_de_novo=False)
   
    annotated_exons_file = filenames.get_exons(config['genome_base_path'], 
        config['genome_name'], note=config.get('orf_note'), 
        is_annotated=True, is_de_novo=False)

    orfs_genomic = filenames.get_orfs(config['genome_base_path'], 
        config['genome_name'], note=config.get('orf_note'))

    exons_file = filenames.get_exons(config['genome_base_path'], 
        config['genome_name'], note=config.get('orf_note'))

   
    # now, check if we have a de novo assembly
    if 'de_novo_gtf' in config:
        get_orfs(config['de_novo_gtf'], args, config, is_annotated=False, 
            is_de_novo=True)

        # we need to concat the ORF and exon files
        de_novo_orfs = filenames.get_orfs(config['genome_base_path'], 
            config['genome_name'], note=config.get('orf_note'), is_annotated=False,
            is_de_novo=True)
       
        de_novo_exons_file = filenames.get_exons(config['genome_base_path'], 
            config['genome_name'], note=config.get('orf_note'), 
            is_annotated=False, is_de_novo=True)

        orfs_files = [annotated_orfs, de_novo_orfs]


        orfs_files_str = ' '.join(orfs_files)
        msg = ("Concatenating files. Output file: {}; Input files: {}".format(
            orfs_genomic, orfs_files_str))
        logger.info(msg)

        if call:
            concatenated_bed = bed_utils.concatenate(orfs_files, sort_bed=True)
            concatenated_bed['orf_num'] = range(len(concatenated_bed))
            fields = bed_utils.bed12_field_names + ['orf_num', 'orf_len', 'orf_type']
            bed_utils.write_bed(concatenated_bed[fields], orfs_genomic)
        else:
            msg = "Skipping concatenation due to --call value"
            logger.info(msg)

        
        exons_files = [annotated_exons_file, de_novo_exons_file]
        
        exons_files_str = ' '.join(exons_files)
        msg = ("Concatenating files. Output file: {}; Input files: {}".format(
            exons_file, exons_files_str))
        logger.info(msg)

        if call:
            concatenated_bed = bed_utils.concatenate(exons_files, sort_bed=True)
            fields = bed_utils.bed6_field_names + ['exon_index', 'transcript_start']
            bed_utils.write_bed(concatenated_bed[fields], exons_file)
        else:
            msg = "Skipping concatenation due to --call value"
            logger.info(msg)

    else:
        # finally, make sure our files are named correctly
        
        if os.path.exists(annotated_orfs):
            utils.create_symlink(annotated_orfs, orfs_genomic, call)

        if os.path.exists(annotated_exons_file):
            utils.create_symlink(annotated_exons_file, exons_file, call)

if __name__ == '__main__':
    main()
