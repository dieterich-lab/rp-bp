#! /usr/bin/env python3

import argparse
import yaml
import logging
import pandas as pd
import os

import numpy as np
import scipy.io

import pbio.ribo.ribo_filenames as ribo_filenames
import pbio.ribo.ribo_utils as ribo_utils


import pbio.utils.bio as bio
import pbio.utils.bam_utils as bam_utils
import pbio.utils.fastx_utils as fastx_utils
import pbio.misc.logging_utils as logging_utils
import pbio.misc.parallel as parallel
import pbio.misc.shell_utils as shell_utils
import pbio.misc.utils as utils
import pbio.misc.pandas_utils as pandas_utils

logger = logging.getLogger(__name__)

default_num_cpus = 2

def get_counts(name_data, config, args):
    name, data = name_data
    msg = "processing {}...".format(name)
    logger.info(msg)

    note = config.get('note', None)
    
    # keep multimappers?
    is_unique = not ('keep_riboseq_multimappers' in config)

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
    unique_bam = ribo_filenames.get_riboseq_bam(
        config['riboseq_data'], name, is_unique=is_unique, note=note)

    # now count the reads of each type
    msg = "{}: collecting read counts".format(name)
    logger.info(msg)
    
    # get the read counts
    msg = "{}: counting reads in raw data".format(name)
    logger.info(msg)
    raw_data_count = fastx_utils.get_read_count(raw_data, is_fasta=False)

    msg = "{}: counting reads without adapters".format(name)
    logger.info(msg)
    without_adapters_count = fastx_utils.get_read_count(without_adapters, is_fasta=False)
    
    msg = "{}: counting reads with rrna".format(name)
    logger.info(msg)
    with_rrna_count = fastx_utils.get_read_count(with_rrna, is_fasta=False)
    
    msg = "{}: counting reads without rrna".format(name)
    logger.info(msg)
    without_rrna_count = fastx_utils.get_read_count(without_rrna, is_fasta=False)
    
    msg = "{}: counting genome-aligned reads".format(name)
    logger.info(msg)
    genome_count = bam_utils.count_aligned_reads(genome_bam)

    msg = "{}: counting uniquely-aligned reads".format(name)
    logger.info(msg)
    unique_count = bam_utils.count_aligned_reads(unique_bam)

    # count reads with correct lengths
    msg = "{}: counting reads with selected lengths".format(name)
    logger.info(msg)

    # now count the unique reads with the appropriate length
    try:
        lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(config, name, 
            is_unique=is_unique)
        lengths_str = ','.join(lengths)
        length_counts = bam_utils.get_length_distribution(unique_bam)
        lengths = set([int(l) for l in lengths])
        m_lengths = length_counts['length'].isin(lengths)
        length_count = np.sum(length_counts.loc[m_lengths, 'count'])

        
        msg = ("{}: found the following periodic lengths: {}. The number of reads "
            "of these lengths: {}".format(name, lengths_str, length_count))
        logger.debug(msg)
    except ValueError as e:
        msg = ("Encountered a problem counting periodic reads. This probably "
            "means no read lengths were periodic. Error message: {}".format(e))
        logger.warning(msg)
        length_count = 0

    ret = {
        'note': name,
        'raw_data_count': raw_data_count,
        'without_adapters_count': without_adapters_count,
        'without_rrna_count': without_rrna_count,
        'genome_count': genome_count,
        'unique_count': unique_count,
        'length_count': length_count
    }

    return pd.Series(ret)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script collects counts of riboseq reads filtered at each step in "
        "the micropeptide prediction pipeline. It mostly parses fastqc results (using the "
        "crimson python package).")
    parser.add_argument('config', help="The yaml config file")
    parser.add_argument('out', help="The output csv file with the counts")
    parser.add_argument('-p', '--num-cpus', help="The number of processors to use", 
        type=int, default=default_num_cpus)
    parser.add_argument('--overwrite', action='store_true')

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    programs = ['samtools']
    shell_utils.check_programs_exist(programs)

    config = yaml.load(open(args.config), Loader=yaml.FullLoader)

    res = parallel.apply_parallel_iter(config['riboseq_samples'].items(), 
                                        args.num_cpus, 
                                        get_counts, config, args)
    res = [r for r in res if r is not None]
    res_df = pd.DataFrame(res)

    pandas_utils.write_df(res_df, args.out, index=False)
    
if __name__ == '__main__':
    main()

