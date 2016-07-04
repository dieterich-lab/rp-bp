#! /usr/bin/env python3

import argparse
import yaml
import logging
import pandas as pd
import os

import numpy as np
import scipy.io

import crimson.fastqc

import riboutils.ribo_filenames as ribo_filenames
import riboutils.ribo_utils as ribo_utils


import misc.bio as bio
import misc.parallel as parallel
import misc.utils as utils

default_tmp = None

default_num_procs = 2

def get_counts(name_data, config, args):
    name, data = name_data
    msg = "processing {}...".format(name)
    logging.info(msg)

    note = config.get('note', None)

    # first, get the ribo_filenames
    raw_data = data
    without_adapters = ribo_filenames.get_without_adapters_fastq(
        config['riboseq_data'], name, note=note)
    with_rrna = ribo_filenames.get_with_rrna_fastq(
        config['riboseq_data'], name, note=note)
    without_rrna = ribo_filenames.get_without_rrna_fastq(
        config['riboseq_data'], name, note=note)
    genome_bam = ribo_filenames.get_riboseq_bam(
        config['riboseq_data'], name, note=note)
    unique_filename = ribo_filenames.get_riboseq_bam(
        config['riboseq_data'], name, is_unique = True, note=note)

    # now, get the fastqc report ribo_filenames
    raw_data_fastqc = ribo_filenames.get_raw_data_fastqc_data(
        config['riboseq_data'], raw_data)
    without_adapters_fastqc = ribo_filenames.get_without_adapters_fastqc_data(
        config['riboseq_data'], name, note=note)
    with_rrna_fastqc = ribo_filenames.get_with_rrna_fastqc_data(
        config['riboseq_data'], name, note=note)
    without_rrna_fastqc = ribo_filenames.get_without_rrna_fastqc_data(
        config['riboseq_data'], name, note=note)

    genome_bam_fastqc = ribo_filenames.get_riboseq_bam_fastqc_data(
        config['riboseq_data'], name, note=note)
    unique_filename_fastqc = ribo_filenames.get_riboseq_bam_fastqc_data(
        config['riboseq_data'], name, is_unique=True, note=note)

    # create the fastqc reports if they do not exist
    raw_data_fastqc_path = ribo_filenames.get_raw_data_fastqc_path(config['riboseq_data'])
    without_adapters_fastqc_path = ribo_filenames.get_without_adapters_fastqc(config['riboseq_data'])
    with_rrna_fastqc_path = ribo_filenames.get_with_rrna_fastqc(config['riboseq_data'])
    without_rrna_fastqc_path = ribo_filenames.get_without_rrna_fastqc(config['riboseq_data'])
    without_rrna_mapping_fastqc_path = ribo_filenames.get_riboseq_bam_fastqc_path(config['riboseq_data'])

    fastqc_tmp_str = ""
    if args.tmp is not None:
        fastqc_tmp_str = "--dir {}".format(args.tmp)

    msg = "Looking for raw data fastqc report: '{}'".format(raw_data_fastqc)
    logging.debug(msg)
    cmd = "fastqc --outdir {} --extract {} {}".format(raw_data_fastqc_path, 
        raw_data, fastqc_tmp_str)
    in_files = [raw_data]
    out_files = [raw_data_fastqc]

    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite)

    cmd = "fastqc --outdir {} --extract {} {}".format(without_adapters_fastqc_path, 
        without_adapters, fastqc_tmp_str)
    in_files = [without_adapters]
    out_files = [without_adapters_fastqc]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite)

    cmd = "fastqc --outdir {} --extract {} {}".format(with_rrna_fastqc_path, 
        with_rrna, fastqc_tmp_str)
    in_files = [with_rrna]
    out_files = [with_rrna_fastqc]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite)

    cmd = "fastqc --outdir {} --extract {} {}".format(without_rrna_fastqc_path, 
        without_rrna, fastqc_tmp_str)
    in_files = [without_rrna]
    out_files = [without_rrna_fastqc]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite)

    cmd = "fastqc --outdir {} --extract {} {}".format(without_rrna_mapping_fastqc_path, 
        genome_bam, fastqc_tmp_str)
    in_files = [genome_bam]
    out_files = [genome_bam_fastqc]

    msg = "genome_bam: '{}'".format(genome_bam)
    logging.debug(msg)
    
    msg = "genome_bam_fastqc: '{}'".format(genome_bam_fastqc)
    logging.debug(msg)

    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite)

    cmd = "fastqc --outdir {} --extract {} {}".format(without_rrna_mapping_fastqc_path, 
        unique_filename, fastqc_tmp_str)
    in_files = [unique_filename]
    out_files = [unique_filename_fastqc]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite)

    # in some cases, fastqc can fail. make sure all of the reports are present
    all_fastqc_reports = [
        raw_data_fastqc,
        without_adapters_fastqc,
        without_rrna_fastqc,
        genome_bam_fastqc,
        unique_filename_fastqc
    ]

    missing_files = [
        f for f in all_fastqc_reports if not os.path.exists(f)
    ]

    if len(missing_files) > 0:
        msg = "The following fastqc reports were not created correctly:\n"
        msg += '\n'.join(missing_files)
        msg += "\nSkipping these counts."
        logging.warning(msg)

        return None


    # and parse the reports
    msg = "{}: parsing FastQC reports".format(name)
    logging.info(msg)
    raw_data_report = crimson.fastqc.parse(raw_data_fastqc)
    without_adapters_report = crimson.fastqc.parse(without_adapters_fastqc)
    without_rrna_report = crimson.fastqc.parse(without_rrna_fastqc)
    genome_bam_report = crimson.fastqc.parse(genome_bam_fastqc)
    unique_filename_report = crimson.fastqc.parse(unique_filename_fastqc)

    msg = "{}: the unique report is: {}".format(name, unique_filename_fastqc)
    logging.debug(msg)

    # get the read counts
    msg = "{}: counting reads in fastqc reports".format(name)
    logging.info(msg)
    
    raw_data_count = sum(l['Count'] 
        for l in raw_data_report['Sequence Length Distribution']['contents'])
    
    without_adapters_count = sum(l['Count'] 
        for l in without_adapters_report['Sequence Length Distribution']['contents'])
    
    without_rrna_count = sum(l['Count'] 
        for l in without_rrna_report['Sequence Length Distribution']['contents'])

    # for the unique counts, we filter out the duplicates, so we can use the fastqc numbers again
    unique_count = sum(l['Count'] 
        for l in unique_filename_report['Sequence Length Distribution']['contents'])

    # fastqc reports duplicated sequences, so for the genome alignments, we have to use samtools to count
    msg = "{}: counting reads with samtools".format(name)
    logging.info(msg)
    genome_count = bio.count_aligned_reads(genome_bam)

    # count reads with correct lengths
    msg = "{}: counting reads with selected lengths".format(name)
    logging.info(msg)

    # now count the unique reads with the appropriate length
    lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(config, name)
    length_counts = bio.get_length_distribution(unique_filename)
    length_count = sum(length_counts[l] for l in length_counts if str(l) in lengths)
    
    lengths_str = ','.join(lengths)
    msg = ("{}: found the following periodic lengths: {}. The number of reads "
        "of these lengths: {}".format(name, lengths_str, length_count))
    logging.debug(msg)

    msg = ("{}: counting filtered reads does not work correctly due to counting "
        "multiple isoforms. It is disabled.".format(name))
    logging.warning(msg)
    transcriptome_count = 0
    #if os.path.exists(signals_filename):
    #    msg = "{}: counting filtered reads aligned to the transcriptome".format(name)
    #    logging.info(msg)
    #    signals = scipy.io.mmread(signals_filename)
    #    transcriptome_count = signals.sum()
    #else:
    #    msg = "{}: could not find the transcript signals: {}. Not counting".format(name, signals_filename)
    #    logging.warning(msg)
    #    transcriptome_count = 0

    ret = {
        'mapping_info_header': 'mapping_info',
        'note': name,
        'raw_data_count': raw_data_count,
        'without_adapters_count': without_adapters_count,
        'without_rrna_count': without_rrna_count,
        'genome_count': genome_count,
        'unique_count': unique_count,
        'length_count': length_count,
        'transcriptome_count': transcriptome_count
    }

    return pd.Series(ret)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script collects counts of riboseq reads filtered at each step in "
        "the micropeptide prediction pipeline. It mostly parses fastqc results (using the "
        "crimson python package).")
    parser.add_argument('config', help="The yaml config file")
    parser.add_argument('out', help="The output csv file with the counts")
    parser.add_argument('-p', '--num-procs', help="The number of processors to use", 
        type=int, default=default_num_procs)
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--tmp', help="Intermediate files (such as fastqc reports when "
        "they are first generated) will be written here", default=default_tmp)
    
    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    programs = ['samtools',
                'fastqc',
                'java']
    utils.check_programs_exist(programs)

    config = yaml.load(open(args.config))

    res = parallel.apply_parallel_iter(config['riboseq_samples'].items(), 
                                        args.num_procs, 
                                        get_counts, config, args)
    res = [r for r in res if r is not None]
    res_df = pd.DataFrame(res)

    utils.write_df(res_df, args.out, index=False)
    
if __name__ == '__main__':
    main()

