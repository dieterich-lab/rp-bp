#! /usr/bin/env python3

import os
import sys
import yaml
import argparse
import logging

import misc.logging_utils as logging_utils
import misc.shell_utils as shell_utils
import misc.slurm as slurm
import misc.utils as utils

import bio_utils.bio as bio
import bio_utils.bed_utils as bed_utils
import bio_utils.star_utils as star_utils

import riboutils.ribo_filenames as filenames

logger = logging.getLogger(__name__)


def get_orfs(gtf, args, config, is_annotated=False, is_de_novo=False):
    """ This helper function processes a GTF file into its ORFs.
    """
    call = not args.do_not_call
    chr_name_file = os.path.join(config['star_index'], 'chrName.txt')
    chr_name_str = "--chr-name-file {}".format(chr_name_file)

    logging_str = logging_utils.get_logging_options_string(args)
    cpus_str = "--num-cpus {}".format(args.num_cpus)
    add_option_str = ''
    if args.add_trx_match:
        add_option_str = '--add-trx-match'

    # extract a BED12 of the annotated ORFs
    transcript_bed = filenames.get_bed(config['genome_base_path'],
                                       config['genome_name'],
                                       is_merged=False,
                                       is_annotated=is_annotated,
                                       is_de_novo=is_de_novo)

    cmd = ("gtf-to-bed12 {} {} {} {} {}".format(gtf,
                                                transcript_bed,
                                                chr_name_str,
                                                cpus_str,
                                                logging_str))
    in_files = [gtf]
    out_files = [transcript_bed]
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files,
                                   overwrite=args.overwrite, call=call)

    exons_file = filenames.get_exons(config['genome_base_path'],
                                     config['genome_name'],
                                     is_annotated=is_annotated,
                                     is_de_novo=is_de_novo)

    cmd = ("split-bed12-blocks {} {} --num-cpus {} {}".format(transcript_bed,
                                                              exons_file,
                                                              args.num_cpus,
                                                              logging_str))
    in_files = [transcript_bed]
    out_files = [exons_file]
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files,
                                   overwrite=args.overwrite, call=call)

    # extract the transcript fasta
    transcript_fasta = filenames.get_transcript_fasta(config['genome_base_path'],
                                                      config['genome_name'],
                                                      is_annotated=is_annotated,
                                                      is_de_novo=is_de_novo)

    cmd = ("extract-bed-sequences {} {} {} {}".format(transcript_bed,
                                                      config['fasta'],
                                                      transcript_fasta,
                                                      logging_str))
    in_files = [transcript_bed, config['fasta']]
    out_files = [transcript_fasta]
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files,
                                   overwrite=args.overwrite, call=call)

    # extract ORFs using genomic coordinates
    orfs_genomic = filenames.get_orfs(config['genome_base_path'],
                                      config['genome_name'],
                                      note=config.get('orf_note'),
                                      is_annotated=is_annotated,
                                      is_de_novo=is_de_novo)
    start_codons_str = utils.get_config_argument(config, 'start_codons')
    stop_codons_str = utils.get_config_argument(config, 'stop_codons')

    cmd = "extract-orf-coordinates {} {} {} {} {} {} {}".format(transcript_bed,
                                                                transcript_fasta,
                                                                orfs_genomic,
                                                                cpus_str,
                                                                start_codons_str,
                                                                stop_codons_str,
                                                                logging_str,
                                                                add_option_str)
    in_files = [transcript_fasta, transcript_bed]
    out_files = [orfs_genomic]
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files,
                                   overwrite=args.overwrite, call=call)

    exons_file = filenames.get_exons(config['genome_base_path'],
                                     config['genome_name'],
                                     note=config.get('orf_note'),
                                     is_annotated=is_annotated,
                                     is_de_novo=is_de_novo,
                                     is_orf=True)

    cmd = ("split-bed12-blocks {} {} --num-cpus {} {}".format(orfs_genomic,
                                                              exons_file,
                                                              args.num_cpus,
                                                              logging_str))
    in_files = [orfs_genomic]
    out_files = [exons_file]
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files,
                                   overwrite=args.overwrite, call=call)

    labeled_orfs = orfs_genomic
    annotated_bed = filenames.get_bed(config['genome_base_path'],
                                      config['genome_name'],
                                      is_merged=False,
                                      is_annotated=True)

    de_novo_str = ""
    if is_de_novo:
        de_novo_str = '--label-prefix "novel_" --filter --nonoverlapping-label "novel" -s'

    cmd = "label-orfs {} {} {} {} {} {} {}".format(annotated_bed,
                                                   orfs_genomic,
                                                   exons_file,
                                                   labeled_orfs,
                                                   cpus_str,
                                                   de_novo_str,
                                                   logging_str)
    in_files = [annotated_bed, orfs_genomic, exons_file]
    #  since we are reusing the name, it will already exist
    out_files = None  # [] # [labeled_orfs]
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files,
                                   overwrite=args.overwrite, call=call)


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='''This script creates all of the files 
                                     necessary for downstream analysis performed with 
                                     the rpbp package.''')

    parser.add_argument('config', help='''The (yaml) config file''')

    parser.add_argument('--overwrite', help='''If this flag is present, existing files
        will be overwritten.''', action='store_true')

    parser.add_argument('--add-trx-match', help='''If this flag is present, an additional
        column is added to the ORFs file containing for each ORF a list of annotated 
        transcripts to which it belongs. This is not used as part of the Rp-Bp pipeline, 
        but may be useful for downstream analysis, this however significantly increase 
        the running time to extract the ORFs.''', action='store_true')

    star_utils.add_star_options(parser)
    slurm.add_sbatch_options(parser)
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    config = yaml.load(open(args.config))
    call = not args.do_not_call

    # check required callable programs, config keys and files
    programs = ['extract-orf-coordinates',
                'label-orfs',
                'bowtie2-build-s',
                'split-bed12-blocks',
                'gtf-to-bed12',
                args.star_executable]
    shell_utils.check_programs_exist(programs)

    required_keys = ['genome_base_path',
                     'genome_name',
                     'gtf',
                     'fasta',
                     'ribosomal_fasta',
                     'ribosomal_index',
                     'star_index']
    utils.check_keys_exist(config, required_keys)

    files = [config['gtf'],
             config['fasta'],
             config['ribosomal_fasta']]
    if 'de_novo_gtf' in config:
        files += [config['de_novo_gtf']]
    utils.check_files_exist(files, source='prepare-rpbp-genome')

    # check if we want to use slurm
    if args.use_slurm:
        cmd = ' '.join(sys.argv)
        slurm.check_sbatch(cmd, args=args)
        return

    # the rRNA index
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
                                                                config['star_index'],
                                                                config['fasta'],
                                                                args.num_cpus,
                                                                mem))

    in_files = [config['fasta']]
    out_files = star_utils.get_star_index_files(config['star_index'])
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files,
                                   overwrite=args.overwrite, call=call)

    # get the ORFs
    get_orfs(config['gtf'], args, config, is_annotated=True, is_de_novo=False)

    # eventually, we will use these names
    annotated_orfs = filenames.get_orfs(config['genome_base_path'],
                                        config['genome_name'], note=config.get('orf_note'), is_annotated=True,
                                        is_de_novo=False)

    annotated_exons_file = filenames.get_exons(config['genome_base_path'],
                                               config['genome_name'], note=config.get('orf_note'),
                                               is_annotated=True, is_de_novo=False, is_orf=True)

    orfs_genomic = filenames.get_orfs(config['genome_base_path'],
                                      config['genome_name'], note=config.get('orf_note'))

    exons_file = filenames.get_exons(config['genome_base_path'],
                                     config['genome_name'], note=config.get('orf_note'), is_orf=True)

    use_gff3_specs = config['gtf'].endswith('gff')
    gtf_file = filenames.get_gtf(config['genome_base_path'],
                                 config['genome_name'], is_gff3=use_gff3_specs, is_star_input=True)

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
                                                 is_annotated=False, is_de_novo=True, is_orf=True)

        orfs_files = [annotated_orfs, de_novo_orfs]

        orfs_files_str = ' '.join(orfs_files)
        msg = ("Concatenating files. Output file: {}; Input files: {}".format(
            orfs_genomic, orfs_files_str))
        logger.info(msg)

        if call:
            concatenated_bed = bed_utils.concatenate(orfs_files, sort_bed=True)
            concatenated_bed['orf_num'] = range(len(concatenated_bed))
            additional_columns = ['orf_num', 'orf_len', 'orf_type']
            if 'assoc_trx' in concatenated_bed.columns:
                additional_columns.extend(['assoc_trx'])
            fields = bed_utils.bed12_field_names + additional_columns
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

        # we also need to concat the annotations to inform STAR
        # there is no particular reason to merge and sort the files, so
        # we just concatenate them...
        if (config['de_novo_gtf'].endswith('gff') == use_gff3_specs):
            cmd = ("awk '!/^#/' {} {} > {}".format(config['gtf'], config['de_novo_gtf'], gtf_file))
            in_files = [config['gtf'], config['de_novo_gtf']]
            out_files = [gtf_file]
            shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files,
                                           overwrite=args.overwrite, call=call)
        else:
            msg = ("Skipping concatenation due to mismatch in format specifications (GTF2/GFF3)"
                   "for reference and do novo annotations. Symlink to reference annotations created.")
            logger.warning(msg)
            if os.path.exists(config['gtf']):
                shell_utils.create_symlink(config['gtf'], gtf_file, call)

    else:
        # finally, make sure our files are named correctly

        if os.path.exists(annotated_orfs):
            shell_utils.create_symlink(annotated_orfs, orfs_genomic, call)

        if os.path.exists(annotated_exons_file):
            shell_utils.create_symlink(annotated_exons_file, exons_file, call)

        if os.path.exists(config['gtf']):
            shell_utils.create_symlink(config['gtf'], gtf_file, call)


if __name__ == '__main__':
    main()
