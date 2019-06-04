#! /usr/bin/env python3

import argparse
import logging
import os
import shlex
import sys

import yaml

import pbio.ribo.ribo_filenames as filenames

import pbio.utils.bio as bio
import pbio.utils.bam_utils as bam_utils
import pbio.utils.fastx_utils as fastx_utils
import pbio.utils.star_utils as star_utils
import pbio.misc.logging_utils as logging_utils
import pbio.misc.shell_utils as shell_utils
import pbio.misc.utils as utils

logger = logging.getLogger(__name__)

default_num_cpus = 1
default_mem = "2G"

# Flexbar arguments
default_max_uncalled = 1
default_pre_trim_left = 0
default_qtrim = "sanger"

flexbar_compression_str = "--zip-output GZ"

# STAR arguments
default_star_executable = "STAR"

default_align_intron_min = 20
default_align_intron_max = 100000
default_out_filter_mismatch_n_max = 1
default_out_filter_mismatch_n_over_l_max = 0.04
default_out_filter_type = "BySJout"
default_out_filter_intron_motifs = "RemoveNoncanonicalUnannotated"
default_out_sam_attributes = ["AS", "NH", "HI", "nM", "MD"]
default_sjdb_overhang = 50

# the Rp-Bp pipeline does not use the transcript alignments, so do not create them
star_out_str = "--outSAMtype BAM SortedByCoordinate"

default_tmp = None


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="")  # filenames.run_riboseq_preprocessing_description)

    parser.add_argument('raw_data', help="The raw data file (fastq[.gz])")
    parser.add_argument('config', help="The (yaml) config file")
    parser.add_argument('name', help="The name for the dataset, used in the created files")

    parser.add_argument('-p', '--num-cpus', help="The number of processors to use",
                        type=int, default=default_num_cpus)
    
    parser.add_argument('--mem', help="The amount of RAM to request",
                        default=default_mem)

    parser.add_argument('--flexbar-options', help="""Optional argument: a space-delimited 
        list of options to pass to flexbar. Each option must be quoted separately as in 
        "--flexbarOption value", using hard, then soft quotes, where "--flexbarOption" 
        is the long parameter name from flexbar and "value" is the value given to this parameter. 
        If specified, flexbar options will override default settings.""", nargs='*', type=str)

    parser.add_argument('-t', '--tmp', help="""The location for temporary files. If not
        specified, program-specific temp locations are used.""", default=default_tmp)

    parser.add_argument('--do-not-call', action='store_true')
    parser.add_argument('--overwrite', help="""If this flag is present, existing files
        will be overwritten.""", action='store_true')

    parser.add_argument('-k', '--keep-intermediate-files', help="""If this flag is given,
        then all intermediate files will be kept; otherwise, they will be
        deleted. This feature is implemented piecemeal. If the --do-not-call flag
        is given, then nothing will be deleted.""", action='store_true')

    star_utils.add_star_options(parser)
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "[create-base-genome-profile]: {}".format(' '.join(sys.argv))
    logger.info(msg)

    config = yaml.load(open(args.config), Loader=yaml.FullLoader)
    call = not args.do_not_call
    keep_delete_files = args.keep_intermediate_files or args.do_not_call

    # check that all of the necessary programs are callable
    programs = [
        'flexbar',
        args.star_executable,
        'samtools',
        'bowtie2',
        'remove-multimapping-reads'
    ]
    shell_utils.check_programs_exist(programs)

    required_keys = [
        'riboseq_data',
        'ribosomal_index',
        'gtf',
        'genome_base_path',
        'genome_name'
    ]
    utils.check_keys_exist(config, required_keys)

    note = config.get('note', None)

    # Step 0: Running flexbar to remove adapter sequences

    raw_data = args.raw_data
    flexbar_target = filenames.get_without_adapters_base(config['riboseq_data'], args.name, note=note)
    without_adapters = filenames.get_without_adapters_fastq(config['riboseq_data'], args.name, note=note)

    adapter_seq_str = utils.get_config_argument(config, 'adapter_sequence', 'adapter-seq')
    adapter_file_str = utils.get_config_argument(config, 'adapter_file', 'adapters')

    max_uncalled_str = utils.get_config_argument(config, 'max_uncalled', default=default_max_uncalled)
    pre_trim_left_str = utils.get_config_argument(config, 'pre_trim_left', default=default_pre_trim_left)
    quality_format_str = utils.get_config_argument(config, "qtrim_format", default=default_qtrim)

    flexbar_option_str = ""
    if args.flexbar_options:
        flexbar_option_str = "{}".format(' '.join(flx_op.strip('"') for flx_op in args.flexbar_options))
    if 'max-uncalled' not in flexbar_option_str:
        flexbar_option_str = "{}".format(' '.join([flexbar_option_str, max_uncalled_str]))
    if 'pre-trim-left' not in flexbar_option_str:
        flexbar_option_str = "{}".format(' '.join([flexbar_option_str, pre_trim_left_str]))
    if 'format' not in flexbar_option_str:
        flexbar_option_str = "{}".format(' '.join([flexbar_option_str, quality_format_str]))

    cmd = "flexbar -r {} -t {} {} {} {} {} -n {}".format(raw_data, flexbar_target, flexbar_compression_str,
                                                         adapter_seq_str, adapter_file_str,
                                                         flexbar_option_str, args.num_cpus)
    in_files = [raw_data]
    out_files = [without_adapters]
    file_checkers = {
        without_adapters: fastx_utils.check_fastq_file
    }
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files,
                                   file_checkers=file_checkers, overwrite=args.overwrite, call=call)

    # Step 1: Running bowtie2 to remove rRNA alignments
    out = utils.abspath("dev", "null")  # we do not care about the alignments
    without_rrna = filenames.get_without_rrna_fastq(config['riboseq_data'], args.name, note=note)
    with_rrna = filenames.get_with_rrna_fastq(config['riboseq_data'], args.name, note=note)

    cmd = "bowtie2 -p {} --very-fast -x {} -U {} -S {} --un-gz {} --al-gz {}".format(
        args.num_cpus, config['ribosomal_index'], without_adapters, out, 
        without_rrna, with_rrna)
    in_files = [without_adapters]
    in_files.extend(bio.get_bowtie2_index_files(config['ribosomal_index']))
    out_files = [without_rrna, with_rrna]
    to_delete = [without_adapters]
    file_checkers = {
        without_rrna: fastx_utils.check_fastq_file
    }
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files,
                                   file_checkers=file_checkers, overwrite=args.overwrite, call=call,
                                   keep_delete_files=keep_delete_files, to_delete=to_delete)

    # Step 2: Running STAR to align rRNA-depleted reads to genome
    star_output_prefix = filenames.get_riboseq_bam_base(config['riboseq_data'], args.name, note=note)
    # transcriptome_bam = "{}{}".format(star_output_prefix, "Aligned.toTranscriptome.out.bam")
    genome_star_bam = "{}{}".format(star_output_prefix, "Aligned.sortedByCoord.out.bam")

    # Get all default options, and add additional options if present.
    # Additional options will override defaults/config options, except for options further below.
    pre_defined_star_options = {}

    align_intron_min_str = utils.get_config_argument(config, 'align_intron_min', 'alignIntronMin',
                                                     default=default_align_intron_min)
    pre_defined_star_options['alignIntronMin'] = align_intron_min_str
    align_intron_max_str = utils.get_config_argument(config, 'align_intron_max', 'alignIntronMax',
                                                     default=default_align_intron_max)
    pre_defined_star_options['alignIntronMax'] = align_intron_max_str
    out_filter_mismatch_n_max_str = utils.get_config_argument(config, 'out_filter_mismatch_n_max',
                                                              'outFilterMismatchNmax',
                                                              default=default_out_filter_mismatch_n_max)
    pre_defined_star_options['outFilterMismatchNmax'] = out_filter_mismatch_n_max_str
    out_filter_mismatch_n_over_l_max_str = utils.get_config_argument(config, 'out_filter_mismatch_n_over_l_max',
                                                                     'outFilterMismatchNoverLmax',
                                                                     default=default_out_filter_mismatch_n_over_l_max)
    pre_defined_star_options['outFilterMismatchNoverLmax'] = out_filter_mismatch_n_over_l_max_str
    out_filter_type_str = utils.get_config_argument(config, 'out_filter_type', 'outFilterType',
                                                    default=default_out_filter_type)
    pre_defined_star_options['outFilterType'] = out_filter_type_str
    out_filter_intron_motifs_str = utils.get_config_argument(config, 'out_filter_intron_motifs',
                                                             'outFilterIntronMotifs',
                                                             default=default_out_filter_intron_motifs)
    pre_defined_star_options['outFilterIntronMotifs'] = out_filter_intron_motifs_str
    out_sam_attributes_str = utils.get_config_argument(config, 'out_sam_attributes', 'outSAMattributes',
                                                       default=default_out_sam_attributes)
    pre_defined_star_options['outSAMattributes'] = out_sam_attributes_str
    sjdb_overhang_str = utils.get_config_argument(config, 'sjdb_overhang', 'sjdbOverhang',
                                                  default=default_sjdb_overhang)
    pre_defined_star_options['sjdbOverhang'] = sjdb_overhang_str
    star_compression_str = "--readFilesCommand {}".format(shlex.quote(args.star_read_files_command))
    pre_defined_star_options['readFilesCommand'] = star_compression_str
    mem_bytes = utils.human2bytes(args.mem)
    star_mem_str = "--limitBAMsortRAM {}".format(mem_bytes)
    pre_defined_star_options['limitBAMsortRAM'] = star_mem_str

    # add additional arguments and/or override defaults
    all_additional_options_str = ""
    if args.star_additional_options:
        all_additional_options_str = "{}".format(' '.join(star_op.strip('"') for star_op
                                                          in args.star_additional_options))
    for star_option_key, star_option_str in pre_defined_star_options.items():
        if star_option_key not in all_additional_options_str:
            all_additional_options_str = "{}".format(' '.join([all_additional_options_str, star_option_str]))

    # tmp directory not handled via star additional options
    star_tmp_str = ""
    if args.tmp is not None:
        star_tmp_name = str(args.name + "_STARtmp")
        star_tmp_dir = star_utils.create_star_tmp(args.tmp, star_tmp_name)
        star_tmp_str = "--outTmpDir {}".format(star_tmp_dir)

    # If GFF3 specs, then we need to inform STAR.
    # Whether we have de novo or not, the format of "config['gtf']" has precedence.
    sjdb_gtf_tag_str = ""
    use_gff3_specs = config['gtf'].endswith('gff')
    gtf_file = filenames.get_gtf(config['genome_base_path'],
                                 config['genome_name'],
                                 is_gff3=use_gff3_specs,
                                 is_star_input=True)
    if use_gff3_specs:
        sjdb_gtf_tag_str = "--sjdbGTFtagExonParentTranscript Parent"

    cmd = ("{} --runThreadN {} --genomeDir {} --sjdbGTFfile {} {} --readFilesIn {} "
        "{} --outFileNamePrefix {} {} {}".format(args.star_executable,
                                                 args.num_cpus,
                                                 config['star_index'],
                                                 gtf_file,
                                                 sjdb_gtf_tag_str,
                                                 without_rrna,
                                                 all_additional_options_str,
                                                 star_output_prefix,
                                                 star_out_str,
                                                 star_tmp_str))
    in_files = [without_rrna]
    in_files.extend(star_utils.get_star_index_files(config['star_index']))
    # out_files = [transcriptome_bam, genome_star_bam]
    to_delete = [without_rrna]
    out_files = [genome_star_bam]
    file_checkers = {
        # transcriptome_bam: bam_utils.check_bam_file,
        genome_star_bam: bam_utils.check_bam_file
    }
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files,
                                   file_checkers=file_checkers, overwrite=args.overwrite,
                                   call=call, keep_delete_files=keep_delete_files, to_delete=to_delete)
    
    # now, we need to symlink the (genome) STAR output to that expected by the rest of the pipeline
    genome_sorted_bam = filenames.get_riboseq_bam(config['riboseq_data'], args.name, note=note)

    if os.path.exists(genome_star_bam):
        shell_utils.create_symlink(genome_star_bam, genome_sorted_bam, call)
    else:
        msg = ("Could not find the STAR genome bam alignment file. Unless "
               "--do-not-call was given, this is a problem.")
        logger.warning(msg)

    # create the bamtools index
    cmd = "samtools index -b {}".format(genome_sorted_bam)
    shell_utils.check_call(cmd, call=call)

    # check if we want to keep multimappers
    if 'keep_riboseq_multimappers' in config:
        return

    # remove multimapping reads from the genome file
    tmp_str = ""
    if args.tmp is not None:
        tmp_str = "--tmp {}".format(args.tmp)

    unique_genome_filename = filenames.get_riboseq_bam(config['riboseq_data'],
                                                       args.name,
                                                       is_unique=True,
                                                       note=note)

    cmd = "remove-multimapping-reads {} {} {}".format(genome_sorted_bam, 
                                                      unique_genome_filename,
                                                      tmp_str)

    in_files = [genome_sorted_bam]
    out_files = [unique_genome_filename]
    to_delete = [genome_star_bam, genome_sorted_bam]
    file_checkers = {
        unique_genome_filename: bam_utils.check_bam_file
    }
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files,
                                   file_checkers=file_checkers, overwrite=args.overwrite,
                                   call=call, keep_delete_files=keep_delete_files, to_delete=to_delete)


if __name__ == '__main__':
    main()
