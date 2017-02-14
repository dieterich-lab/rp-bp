#! /usr/bin/env python3

import argparse
import shlex
import string
import sys
import yaml

import misc.slurm as slurm
import misc.utils as utils

import riboutils.ribo_filenames as filenames
import riboutils.ribo_utils as ribo_utils

import logging
import misc.logging_utils as logging_utils
logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Extract the ORF profiles for each specified read length "
        "and offset independently. One sparse matrix file will be created for "
        "each read length. It then collects the values into a sparse tensor.")

    parser.add_argument('config', help="The (json) config file")
    parser.add_argument('name', help="The name for the dataset, used in the "
        "created files")
    
    parser.add_argument('out', help="The (mtx.gz) output file containing the "
        "ORF profiles and read lengths")

    parser.add_argument('-c', '--is-condition', help="If this flag is present, "
        "then \"name\" will be taken to be a condition name. The profiles for "
        "all relevant replicates of the condition will be created.", 
        action='store_true')

    slurm.add_sbatch_options(parser)
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    logging_str = logging_utils.get_logging_options_string(args)
    cpus_str = "--num-cpus {}".format(args.num_cpus)

    msg = "[create-read-length-orf-profiles]: {}".format(' '.join(sys.argv))
    logger.info(msg)

    msg = "Reading config file"
    logger.info(msg)
    config = yaml.load(open(args.config))
 
    # pull out what we need from the config file
    is_unique = not ('keep_riboseq_multimappers' in config)    
    seqname_str = utils.get_config_argument(config, 'seqname_prefix')
    note = config.get('note', None)
    orf_note = config.get('orf_note', None)

    
    orfs = filenames.get_orfs(
        config['genome_base_path'], 
        config['genome_name'], 
        note=orf_note
    )

    exons = filenames.get_exons(
        config['genome_base_path'], 
        config['genome_name'],
        note=orf_note
    )
    
    # make sure the necessary files exist
    required_files = [orfs, exons]
    msg = "[create-read-length-orf-profiles]: Some input files were missing: "
    utils.check_files_exist(required_files, msg=msg)

    # check which samples to process
    names = [args.name]
    if args.is_condition:
        riboseq_replicates = ribo_utils.get_riboseq_replicates(config)
        names = [n for n in riboseq_replicates[args.name]]

    job_ids = []
    for name in names:

        msg = "Processing sample: {}".format(name)
        logger.info(msg)
        
        # now the relevant files
        bam = filenames.get_riboseq_bam(
            config['riboseq_data'], 
            name, 
            is_unique=is_unique, 
            note=note
        )

        # make sure the necessary files exist
        required_files = [bam]
        msg = "[create-read-length-orf-profiles]: Some input files were missing: "
        utils.check_files_exist(required_files, msg=msg)

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
            lengths_str = "--lengths {}".format(length)
            offsets_str = "--offsets {}".format(offset)

            
            mtx = filenames.get_riboseq_profiles(
                config['riboseq_data'], 
                name, 
                length=[length], 
                offset=[offset],
                is_unique=is_unique, 
                note=note
            )

            cmd = "extract-orf-profiles {} {} {} {} {} {} {} {}".format(
                bam,
                orfs,
                exons,
                mtx,
                lengths_str,
                offsets_str,
                seqname_str,
                cpus_str,
                logging_str
            )
            
            job_id = slurm.check_sbatch(cmd, args=args)

            job_ids.append(job_id)

    # now, collect them into a single file
    offsets_str = ' '.join(str(o) for o in offsets)
    lengths_str = ' '.join(str(l) for l in lengths)

    offsets_str = "--offsets {}".format(offsets_str)
    lengths_str = "--lengths {}".format(lengths_str)

    is_condition_str = ""
    if args.is_condition:
        is_condition_str = "--is-condition"

    cmd = "collect-read-length-orf-profiles {} {} {} {} {}".format(
        args.config,
        args.name,
        args.out,
        is_condition_str,
        logging_str
    )

    slurm.check_sbatch(cmd, args=args, dependencies=job_ids)

if __name__ == '__main__':
    main()
