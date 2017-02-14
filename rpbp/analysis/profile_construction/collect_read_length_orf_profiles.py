#! /usr/bin/env python3

import argparse
import gzip
import string
import scipy.io
import yaml

import misc.utils as utils

import riboutils.ribo_filenames as filenames
import riboutils.ribo_utils as ribo_utils

import logging
import misc.logging_utils as logging_utils
logger = logging.getLogger(__name__)

default_lengths = []
default_offsets = []

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Collect the individual read length ORF profiles created "
        "by create-read-length-orf-profiles into a single sparse \"tensor\". "
        "N.B. This script is called by create-read-length-orf-profiles, and "
        "it is unlikely that it should ever be called independently.")
    
    parser.add_argument('config', help="The (json) config file")
    parser.add_argument('name', help="The name for the dataset, used in the "
        "created files")
    
    parser.add_argument('out', help="The (mtx.gz) output file containing the "
        "ORF profiles and read lengths")

    parser.add_argument('-c', '--is-condition', help="If this flag is present, "
        "then \"name\" will be taken to be a condition name. The profiles for "
        "all relevant replicates of the condition will be created.", 
        action='store_true')

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)
    
    msg = "Reading config file"
    logger.info(msg)
    config = yaml.load(open(args.config))

    # pull out what we need from the config file
    is_unique = not ('keep_riboseq_multimappers' in config)    
    seqname_prefix_str = utils.get_config_argument(config, 'seqname_prefix')
    note = config.get('note', None)

    names = [args.name]
    if args.is_condition:
        riboseq_replicates = ribo_utils.get_riboseq_replicates(config)
        names = [n for n in riboseq_replicates[args.name]]

    # keep a map from the lengths to the combined profiles
    length_profile_map = {}

    for name in names:

        msg = "Processing sample: {}".format(name)
        logger.info(msg)
        
        # get the lengths and offsets which meet the required criteria from the config file
        lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(
            config, 
            name, 
            is_unique=is_unique
        )

        lengths = [28, 29]
        offsets = [12, 12]

        if len(lengths) == 0:
            msg = ("No periodic read lengths and offsets were found. Try relaxing "
                "min_metagene_profile_count, min_metagene_bf_mean, "
                "max_metagene_bf_var, and/or min_metagene_bf_likelihood. Qutting.")
            logger.critical(msg)
            return

        for length, offset in zip(lengths, offsets):
                        
            mtx = filenames.get_riboseq_profiles(
                config['riboseq_data'], 
                name, 
                length=[length], 
                offset=[offset], 
                is_unique=is_unique, 
                note=note
            )

            mtx = scipy.io.mmread(mtx).tocsr()

            prior_mtx = length_profile_map.get(length, None)

            if prior_mtx is None:
                length_profile_map[length] = mtx
            else:
                length_profile_map[length] = prior_mtx + mtx

    with gzip.open(args.out, 'wb') as target_gz:

        for length, mtx in length_profile_map.items():
            mtx = mtx.tocoo()
                        
            msg = "Writing ORF profiles. length: {}.".format(length)
            logger.info(msg)
            
            for row, col, val in zip(mtx.row, mtx.col, mtx.data):
                s = "{} {} {} {}\n".format(length, row, col, val)
                target_gz.write(s.encode())

if __name__ == '__main__':
    main()
