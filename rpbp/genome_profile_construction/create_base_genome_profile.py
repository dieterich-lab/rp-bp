#! /usr/bin/env python3

import argparse
import logging
import os

import yaml

import rpbp.filenames as filenames
import misc.utils as utils

default_num_procs = 2

default_quality_format = 'sanger'
default_max_uncalled = 1
default_pre_trim_left = 0

# STAR arguments
default_star_executable = "STAR"

default_align_intron_min = 20
default_align_intron_max = 100000
default_out_filter_mismatch_n_max = 1
default_out_filter_mismatch_n_over_l_max = 0.04
default_out_filter_type = "BySJout"
default_out_filter_intron_motifs = "RemoveNoncanonicalUnannotated"
default_out_sam_attributes = "AS NH HI nM MD"

flexbar_compression_str = "--zip-output GZ"
quant_mode_str = '--quantMode TranscriptomeSAM'
star_compression_str = "--readFilesCommand zcat"
star_out_str = "--outSAMtype BAM SortedByCoordinate"

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=filenames.run_riboseq_preprocessing_description)

    parser.add_argument('raw_data', help="The raw data file (fastq[.gz])")
    parser.add_argument('config', help="The (yaml) config file")
    parser.add_argument('name', help="The name for the dataset, used in the created files")

    parser.add_argument('-p', '--num-procs', help="The number of processors to use",
        type=int, default=default_num_procs)
        
    parser.add_argument('--star-executable', help="The name of the STAR executable",
        default=default_star_executable)

    parser.add_argument('--do-not-call', action='store_true')
    parser.add_argument('--overwrite', help="If this flag is present, existing files "
        "will be overwritten.", action='store_true')
    
    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    config = yaml.load(open(args.config))
    call = not args.do_not_call

    # check that all of the necessary programs are callable
    programs =  [   'flexbar',
                    args.star_executable,
                    'samtools',
                    'bowtie2',
                    'remove-multimapping-reads'
                ]
    utils.check_programs_exist(programs)

    required_keys = [   'riboseq_data',
                        'ribosomal_index',
                        'gtf',
                        'star_index'
                    ]
    utils.check_keys_exist(config, required_keys)


    # Step 0: Running flexbar to remove adapter sequences

    raw_data = args.raw_data
    flexbar_target = filenames.get_without_adapters_base(config['riboseq_data'], args.name)
    without_adapters = filenames.get_without_adapters_fastq(config['riboseq_data'], args.name)

    adapter_seq_str = utils.get_config_argument(config, 'adapter_sequence', 'adapter-seq')
    adapter_file_str = utils.get_config_argument(config, 'adapter_file', 'adapters')

    quality_format_str = utils.get_config_argument(config, 'quality_format', 'format', 
        default=default_quality_format)
    max_uncalled_str = utils.get_config_argument(config, 'max_uncalled', default=default_max_uncalled)
    pre_trim_left_str = utils.get_config_argument(config, 'pre_trim_left', default=default_pre_trim_left)

    cmd = "flexbar {} {} {} {} -n {} {} -r {} -t {} {}".format(quality_format_str, 
        max_uncalled_str, adapter_seq_str, adapter_file_str, args.num_procs, flexbar_compression_str, 
        raw_data, flexbar_target, pre_trim_left_str)
    in_files = [raw_data]
    out_files = [without_adapters]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)

    # Step 1: Running bowtie2 to remove rRNA alignments
    out = utils.abspath("dev","null") # we do not care about the alignments
    without_rrna = filenames.get_without_rrna_fastq(config['riboseq_data'], args.name)
    with_rrna = filenames.get_with_rrna_fastq(config['riboseq_data'], args.name)

    cmd = "bowtie2 -p {} --very-fast -x {} -U {} -S {} --un-gz {} --al-gz {}".format(
        args.num_procs, config['ribosomal_index'], without_adapters, out, 
        without_rrna, with_rrna)
    in_files = [without_adapters]
    out_files = [without_rrna, with_rrna]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)

    # Step 2: Running STAR to align rRNA-depleted reads to genome
    star_output_prefix = filenames.get_riboseq_sam_base(config['riboseq_data'], args.name)
    transcriptome_bam = "{}{}".format(star_output_prefix, "Aligned.toTranscriptome.out.bam")
    genome_star_bam = "{}{}".format(star_output_prefix, "Aligned.sortedByCoord.out.bam")

    align_intron_min_str = utils.get_config_argument(config, 'align_intron_min', 
        'alignIntronMin', default=default_align_intron_min)
    align_intron_max_str = utils.get_config_argument(config, 'align_intron_max', 
        'alignIntronMax', default=default_align_intron_max)
    out_filter_mismatch_n_max_str = utils.get_config_argument(config, 'out_filter_mismatch_n_max', 
        'outFilterMismatchNmax', default=default_out_filter_mismatch_n_max)
    out_filter_mismatch_n_over_l_max_str = utils.get_config_argument(config, 'out_filter_mismatch_n_over_l_max',
        'outFilterMismatchNoverLmax', default=default_out_filter_mismatch_n_over_l_max)
    out_filter_type_str = utils.get_config_argument(config, 'out_filter_type', 
        'outFilterType', default=default_out_filter_type)
    out_filter_intron_motifs_str = utils.get_config_argument(config, 'out_filter_intron_motifs', 
        'outFilterIntronMotifs', default=default_out_filter_intron_motifs)
    out_sam_attributes_str = utils.get_config_argument(config, 'out_sam_attributes', 
        'outSAMattributes', default=default_out_sam_attributes)

    cmd = ("{} --runThreadN {} {} --genomeDir {} --sjdbGTFfile {} --readFilesIn {} "
        "{} {} {} {} {} {} {} {} --outFileNamePrefix {} {}".format(args.star_executable,
        args.num_procs, star_compression_str, config['star_index'], config['gtf'], without_rrna, 
        align_intron_min_str, align_intron_max_str, out_filter_mismatch_n_max_str, 
        out_filter_type_str, out_filter_intron_motifs_str, quant_mode_str,
        out_filter_mismatch_n_over_l_max_str, out_sam_attributes_str, star_output_prefix,
        star_out_str))
    in_files = [without_rrna]
    out_files = [transcriptome_bam, genome_star_bam]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)

    # now, we need to rename the STAR output to that expected by the rest of the pipeline
    genome_sorted_bam = filenames.get_riboseq_bam(config['riboseq_data'], args.name)
    os.replace(genome_star_bam, genome_sorted_bam)
     
    # remove the multimapping riboseq reads
    unique_filename = filenames.get_riboseq_bam(config['riboseq_data'], args.name, is_unique=True)
    cmd = "remove-multimapping-reads {} {}".format(genome_sorted_bam, unique_filename)
    in_files = [genome_sorted_bam]
    out_files = [unique_filename]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)

    # Step 3: Running samtools to get the sorted bam files from the STAR transcriptome alignments
    transcriptome_sorted_bam = filenames.get_riboseq_bam(
        config['riboseq_data'], args.name, is_transcriptome=True)
    
    cmd = "samtools view -bS -@{} {} | samtools sort - -@{} -o {}".format(args.num_procs, 
        transcriptome_bam, args.num_procs, transcriptome_sorted_bam)
    in_files = [transcriptome_bam]
    out_files = [transcriptome_sorted_bam]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)

    # Step 4: Creating bam index file for transcriptome alignments
    transcriptome_sorted_bai = transcriptome_sorted_bam + ".bai"
    cmd = "samtools index -b {} {}".format(transcriptome_sorted_bam, transcriptome_sorted_bai)
    in_files = [transcriptome_sorted_bam]
    out_files = [transcriptome_sorted_bai]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call)
   
if __name__ == '__main__':
    main()
