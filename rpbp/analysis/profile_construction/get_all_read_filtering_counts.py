#! /usr/bin/env python3

import argparse
import yaml
import logging
import pandas as pd
import os
import scipy.io

import crimson.fastqc

import riboseq_utils.riboseq_filenames as filenames

import misc.bio as bio
import misc.parallel as parallel
import misc.utils as utils

default_num_procs = 2

def get_counts(name, config, overwrite):
    msg = "processing {}...".format(name)
    logging.info(msg)

    # first, get the filenames
    raw_data = filenames.get_riboseq_fastq(config['riboseq_data'], name)
    without_adapters = filenames.get_without_adapters_fastq(config['riboseq_data'], name)
    with_rrna = filenames.get_with_rrna_fastq(config['riboseq_data'], name)
    without_rrna = filenames.get_without_rrna_fastq(config['riboseq_data'], name)
    genome_bam = filenames.get_riboseq_bam(config['riboseq_data'], name)
    unique_filename = filenames.get_riboseq_bam(config['riboseq_data'], name, is_unique = True)

    # now, get the fastqc report filenames
    raw_data_fastqc = filenames.get_raw_data_fastqc_data(config['riboseq_data'], name)
    without_adapters_fastqc = filenames.get_without_adapters_fastqc_data(config['riboseq_data'], name)
    with_rrna_fastqc = filenames.get_with_rrna_fastqc_data(config['riboseq_data'], name)
    without_rrna_fastqc = filenames.get_without_rrna_fastqc_data(config['riboseq_data'], name)

    genome_bam_fastqc = filenames.get_riboseq_bam_fastqc_data(config['riboseq_data'], name)
    unique_filename_fastqc = filenames.get_riboseq_bam_fastqc_data(config['riboseq_data'], name)

    # create the fastqc reports if they do not exist
    raw_data_fastqc_path = filenames.get_raw_data_fastqc_path(config['riboseq_data'])
    without_adapters_fastqc_path = filenames.get_without_adapters_fastqc(config['riboseq_data'])
    with_rrna_fastqc_path = filenames.get_with_rrna_fastqc(config['riboseq_data'])
    without_rrna_fastqc_path = filenames.get_without_rrna_fastqc(config['riboseq_data'])
    without_rrna_mapping_fastqc_path = filenames.get_riboseq_bam_fastqc_path(config['riboseq_data'])

    cmd = "fastqc --outdir {} --extract {}".format(raw_data_fastqc_path, raw_data)
    in_files = [raw_data]
    out_files = [raw_data_fastqc]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=overwrite)

    cmd = "fastqc --outdir {} --extract {}".format(without_adapters_fastqc_path, without_adapters)
    in_files = [without_adapters]
    out_files = [without_adapters_fastqc]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=overwrite)

    cmd = "fastqc --outdir {} --extract {}".format(with_rrna_fastqc_path, with_rrna)
    in_files = [with_rrna]
    out_files = [with_rrna_fastqc]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=overwrite)

    cmd = "fastqc --outdir {} --extract {}".format(without_rrna_fastqc_path, without_rrna)
    in_files = [without_rrna]
    out_files = [without_rrna_fastqc]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=overwrite)

    cmd = "fastqc --outdir {} --extract {}".format(without_rrna_mapping_fastqc_path, genome_bam)
    in_files = [genome_bam]
    out_files = [genome_bam_fastqc]

    msg = "genome_bam: '{}'".format(genome_bam)
    logging.debug(msg)
    
    msg = "genome_bam_fastqc: '{}'".format(genome_bam_fastqc)
    logging.debug(msg)

    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=overwrite)

    cmd = "fastqc --outdir {} --extract {}".format(without_rrna_mapping_fastqc_path, unique_filename)
    in_files = [unique_filename]
    out_files = [unique_filename_fastqc]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=overwrite)

    # and parse the reports
    msg = "{}: parsing FastQC reports".format(name)
    logging.info(msg)
    raw_data_report = crimson.fastqc.parse(raw_data_fastqc)
    without_adapters_report = crimson.fastqc.parse(without_adapters_fastqc)
    without_rrna_report = crimson.fastqc.parse(without_rrna_fastqc)
    genome_bam_report = crimson.fastqc.parse(genome_bam_fastqc)
    unique_filename_report = crimson.fastqc.parse(unique_filename_fastqc)

    # get the read counts
    msg = "{}: counting reads in fastqc reports".format(name)
    logging.info(msg)
    raw_data_count = sum(l['Count'] for l in raw_data_report['Sequence Length Distribution']['contents'])
    without_adapters_count = sum(l['Count'] for l in without_adapters_report['Sequence Length Distribution']['contents'])
    without_rrna_count = sum(l['Count'] for l in without_rrna_report['Sequence Length Distribution']['contents'])
    # for the unique counts, we filter out the duplicates, so we can use the fastqc numbers again
    unique_count = sum(l['Count'] for l in unique_filename_report['Sequence Length Distribution']['contents'])

    # fastqc reports duplicated sequences, so for the genome alignments, we have to use samtools to count
    msg = "{}: counting reads with samtools".format(name)
    logging.info(msg)
    genome_count = bio.count_aligned_reads(genome_bam)

    # count reads with correct lengths
    msg = "{}: counting reads with selected lengths".format(name)
    logging.info(msg)
    best_periodicity_and_offsets = filenames.get_riboseq_best_periodicity_and_offsets(config['riboseq_data'], 
        name, is_unique=True)
    best_info_df = pd.read_csv(best_periodicity_and_offsets)
    mask_largest_count = best_info_df['largest_count'] > config['min_periodicity_count']
    mask_largest_count_bf = best_info_df['largest_count_bf'] > config['min_periodicity_bf']
    filtered_po = best_info_df[mask_largest_count & mask_largest_count_bf]

    offsets = filtered_po['largest_count_offset']
    lengths = filtered_po['length']
    
    # offsets must be positive
    offsets_l = [str(-1*int(o)) for o in offsets]
    lengths_l = [str(int(l)) for l in lengths]
    
    # now count the unique reads with the appropriate length
    length_count = sum(l['Count'] for l in unique_filename_report['Sequence Length Distribution']['contents'] if str(l['Length']) in lengths_l)

    msg = "{}: counting filtered reads does not work correctly due to counting multiple isoforms. It is disabled.".format(name)
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
    
    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    programs = ['samtools',
                'fastqc',
                'java']
    utils.check_programs_exist(programs)

    config = yaml.load(open(args.config))

    res = parallel.apply_parallel_iter(config['sample_names'], args.num_procs, get_counts, config, args.overwrite)
    res_df = pd.DataFrame(res)
    #res_df = pd.concat(res)

    utils.write_df(res_df, args.out, index=False)
    
if __name__ == '__main__':
    main()

