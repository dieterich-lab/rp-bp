#! /usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import logging
import sys

import tqdm
import pysam

import misc.bio as bio
import misc.parallel as parallel
import misc.utils as utils


default_num_cpus = 1
default_lengths = []
default_num_alignments = 0

default_start_upstream = 50
default_start_downstream = 20
default_end_upstream = 50
default_end_downstream = 20

default_seqids_to_keep = []

def find_forward_first_matching_read_positions(cds_region, reads):
    upstream = cds_region['start_upstream']
    downstream= cds_region['start_downstream']
    
    mask_start = (reads['start'] >= upstream) & (reads['start'] < downstream)
    matching_reads = reads[mask_start]
    relative_positions = np.array(matching_reads['start'] - upstream, dtype=int)
    return relative_positions

def find_forward_last_matching_read_positions(cds_region, reads):
    upstream = cds_region['end_upstream']
    downstream= cds_region['end_downstream']
    
    mask_start = (reads['start'] >= upstream) & (reads['start'] < downstream)
    matching_reads = reads[mask_start]
    relative_positions = np.array(matching_reads['start'] - upstream, dtype=int)
    return relative_positions

def find_reverse_first_matching_read_positions(cds_region, reads):
    upstream = cds_region['start_upstream']
    downstream= cds_region['start_downstream']
    
    mask_start = (reads['start'] <= upstream) & (reads['start'] >= downstream)
    matching_reads = reads[mask_start]
    relative_positions = np.array(matching_reads['start'] - downstream, dtype=int)
    return relative_positions

def find_reverse_last_matching_read_positions(cds_region, reads):
    upstream = cds_region['end_upstream']
    downstream= cds_region['end_downstream']
    
    mask_start = (reads['start'] <= upstream) & (reads['start'] >= downstream)
    matching_reads = reads[mask_start]
    relative_positions = np.array(matching_reads['start'] - downstream, dtype=int)
    return relative_positions

def get_metagene_profile(length_alignment_df, args):
    length, alignment_df = length_alignment_df
    
    logging.debug("Getting annotated start and end regions...")
    # pull out the canonical CDS regions
    orfs = bio.read_bed(args.orfs)

    m_canonical = orfs['orf_type'] == 'canonical'
    m_forward = orfs['strand'] == '+'
    m_reverse = ~m_forward
    
    #canonical_forward_df = orfs.loc[m_canonical & m_forward]
    #canonical_reverse_df = orfs.loc[m_canonical & m_reverse]
    m_canonical_forward = m_canonical & m_forward
    m_canonical_reverse = m_canonical & m_reverse

    msg = "Found {} forward canonical and {} reverse canonical ORFs".format(sum(m_canonical_forward), sum(m_canonical_reverse))
    logging.debug(msg)

    # start and end already contain the boundaries of the ORF

    # TODO: these calls raise a pandas SettingWithCopyWarning. It does not look
    # very good when executing the script.

    # set the ranges we want to find
    #canonical_forward_df['start_upstream'] = canonical_forward_df['start'] - args.start_upstream
    #canonical_forward_df['start_downstream'] = canonical_forward_df['start'] + args.start_downstream
    orfs.loc[m_canonical_forward, 'start_upstream'] = orfs.loc[m_canonical_forward, 'start'] - args.start_upstream
    orfs.loc[m_canonical_forward, 'start_downstream'] = orfs.loc[m_canonical_forward, 'start'] + args.start_downstream

    #canonical_forward_df['end_upstream'] = canonical_forward_df['end'] - args.end_upstream
    #canonical_forward_df['end_downstream'] = canonical_forward_df['end'] + args.end_downstream
    orfs.loc[m_canonical_forward, 'end_upstream'] = orfs.loc[m_canonical_forward, 'end'] - args.end_upstream
    orfs.loc[m_canonical_forward, 'end_downstream'] = orfs.loc[m_canonical_forward, 'end'] + args.end_downstream


    # WE ARE SWITCHING THE ORDER OF THE COORDINATES FOR THE REVERSE STRAND
    # AFTER THIS, start_upstream, etc., WILL HAVE THE SAME SEMANTICS AS FOR THE FORWARD STRAND
    
    #canonical_reverse_df['start_upstream'] = canonical_reverse_df['end'] + args.start_upstream
    #canonical_reverse_df['start_downstream'] = canonical_reverse_df['end'] - args.start_downstream
    orfs.loc[m_canonical_reverse, 'start_upstream'] = orfs.loc[m_canonical_reverse, 'end'] + args.start_upstream
    orfs.loc[m_canonical_reverse, 'start_downstream'] = orfs.loc[m_canonical_reverse, 'end'] - args.start_downstream

    #canonical_reverse_df['end_upstream'] = canonical_reverse_df['start'] + args.end_upstream
    #canonical_reverse_df['end_downstream'] = canonical_reverse_df['start'] - args.end_downstream
    orfs.loc[m_canonical_reverse, 'end_upstream'] = orfs.loc[m_canonical_reverse, 'start'] + args.end_upstream
    orfs.loc[m_canonical_reverse, 'end_downstream'] = orfs.loc[m_canonical_reverse, 'start'] - args.end_downstream
    
    forward_first_cds_count = np.zeros(args.start_upstream + args.start_downstream + 1)
    forward_last_cds_count = np.zeros(args.end_upstream + args.end_downstream + 1)

    reverse_first_cds_count = np.zeros(args.start_upstream + args.start_downstream + 1)
    reverse_last_cds_count = np.zeros(args.end_upstream + args.end_downstream + 1)
    
    if len(args.seqids_to_keep) > 0:
        seqnames = args.seqids_to_keep
    else:
        #seqnames = canonical_forward_df['seqname'].unique()
        seqnames = orfs.loc[m_canonical_forward, 'seqname'].unique()

    for seqname in seqnames:
        msg = "seqname: {}".format(seqname)
        logging.debug(msg)
        
        mask_reads_seq = alignment_df['seqname'] == seqname
        reads_region = alignment_df[mask_reads_seq]
        
        m_orf_seq = orfs['seqname'] == seqname
        
        # first, handle the forward, first CDS reads
        logging.debug("ff...")
        #mask_ff_seq = orfs['seqname'] == seqname
        #forward_seq_df = canonical_forward_df[mask_ff_seq]
        forward_seq_df = orfs[m_orf_seq & m_canonical_forward]

        res = parallel.apply_df_simple(forward_seq_df, find_forward_first_matching_read_positions, reads_region)
        if len(res) > 0:
            res = np.concatenate(res)
        for i in res:
            forward_first_cds_count[i] += 1

        msg = "Found {} reads in region".format(len(res))
        logging.debug(msg)
            
        # the forward, last CDS reads
        logging.debug("fl...")
        res = parallel.apply_df_simple(forward_seq_df, find_forward_last_matching_read_positions, reads_region)
        if len(res) > 0:
            res = np.concatenate(res)
        for i in res:
            forward_last_cds_count[i] += 1

        msg = "Found {} reads in region".format(len(res))
        logging.debug(msg)
            
        # reverse, first CDS reads
        logging.debug("rf...")
        #m_rf_seq = canonical_reverse_df['seqname'] == seqname
        #reverse_seq_df = canonical_reverse_df[mask_rf_seq]
        reverse_seq_df = orfs[m_orf_seq & m_canonical_reverse]
        res = parallel.apply_df_simple(reverse_seq_df, find_reverse_first_matching_read_positions, reads_region)
        if len(res) > 0:
            res = np.concatenate(res)
        for i in res:
            reverse_first_cds_count[i] += 1

        msg = "Found {} reads in region".format(len(res))
        logging.debug(msg)
            
        # reverse, last CDS reads
        logging.debug("rl...")
        res = parallel.apply_df_simple(reverse_seq_df, find_reverse_last_matching_read_positions, reads_region)
        if len(res) > 0:
            res = np.concatenate(res)
        for i in res:
            reverse_last_cds_count[i] += 1

        msg = "Found {} reads in region".format(len(res))
        logging.debug(msg)

    # we need to reverse the "reverse" counts so they match the forward strand counts
    reverse_first_cds_count_r = reverse_first_cds_count[::-1]
    reverse_last_cds_count_r = reverse_last_cds_count[::-1]

    start_counts = reverse_first_cds_count_r + forward_first_cds_count
    end_counts = reverse_last_cds_count_r + forward_last_cds_count

    start_positions = range(-1*args.start_upstream, args.start_downstream + 1)
    end_positions = range(-1*args.end_upstream, args.end_downstream + 1)

    logging.debug("len(start_positions): {}".format(len(start_positions)))
    logging.debug("len(start_counts): {}".format(len(start_counts)))

    # now, we need to dump this signal out
    start_df = pd.DataFrame()
    start_df['count'] = start_counts
    start_df['position'] = start_positions
    start_df['type'] = 'start'

    end_df = pd.DataFrame()
    end_df['count'] = end_counts
    end_df['position'] = end_positions
    end_df['type'] = 'end'

    length_profile_df = pd.concat([start_df, end_df])
    length_profile_df['length'] = length
    
    return length_profile_df


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script extracts the metagene profile from reads in a BAM "
        "file, possibly filtering by length. It attempts to vectorize as many of the "
        "counting operations as possible.")
    parser.add_argument('bam', help="The bam file")
    parser.add_argument('orfs', help="The annotated orfs (bed) file")
    parser.add_argument('out', help="The (output) csv.gz counts file")

    parser.add_argument('-p', '--num-cpus', help="The number of processors to use",
        type=int, default=default_num_cpus)

    parser.add_argument('--is-sam', help="If this flag is present, the alignment file will "
        "be parsed as SAM rather than BAM", action='store_true')

    parser.add_argument('--lengths', help="If specified, then metagene profiles will be "
        "created for reads of each length. Otherwise, profiles will be created for each "
        "read length present in the bam file.", type=int, nargs='*', default=default_lengths)

    parser.add_argument('--seqids-to-keep', help="If any seqids are given here, then "
        "only those will be retained.", nargs='+', default=default_seqids_to_keep)

    parser.add_argument('--num-alignments', help="If this value is >0, then only the "
        "first k alignments in the file will be analyzes. This is implemented using "
        "pysam.AlignmentFile.head, so its semantics determine the analyzed alignments.",
        type=int, default=default_num_alignments)

    parser.add_argument('--start-upstream', type=int, default=default_start_upstream, 
        help="The number of bases upstream of the translation initiation site to begin "
        "constructing the metagene profile.")
    parser.add_argument('--start-downstream', type=int, default=default_start_downstream,
        help="The number of bases downstream of the translation initiation site to end "
        "the metagene profile.")
    parser.add_argument('--end-upstream', type=int, default=default_end_upstream,
        help="The number of bases upstream of the translation termination site to begin "
        "constructing the metagene profile.")
    parser.add_argument('--end-downstream', type=int, default=default_end_downstream,
        help="The number of bases downstream of the translation termination site to end "
        "the metagene profile.")

    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    if args.is_sam:
        bam = pysam.AlignmentFile(args.bam, 'r')
    else:
        bam = pysam.AlignmentFile(args.bam, 'rb')

    logging.info("Processing alignments...")

    if args.num_alignments != default_num_alignments:
        alignments = bam.head(args.num_alignments)
        num_alignments = args.num_alignments
    else:
        alignments = bam.fetch()
        num_alignments = bam.count()
    
    lengths = np.zeros(num_alignments)
    starts = np.zeros(num_alignments)
    seqs = [None] * num_alignments

    for i, a in enumerate(tqdm.tqdm(alignments, leave=True, file=sys.stdout, total=num_alignments)):
        starts[i] = a.reference_start
        if a.is_reverse:
            starts[i] = a.reference_end
            
        lengths[i] = a.qlen
        seqs[i] = a.reference_name

    alignment_df = pd.DataFrame()
    alignment_df['start'] = starts
    alignment_df['length'] = lengths
    alignment_df['seqname'] = seqs
    
    msg = "Reads remaining after filtering: {}".format(len(alignment_df))
    logging.info(msg)

    if len(args.lengths) == 0:
        args.lengths = list(alignment_df['length'].unique())

    length_str = ','.join(str(int(l)) for l in args.lengths)
    msg = "Profiles will be created for lengths: {}".format(length_str)
    logging.info(msg)


    # now, split out the individual length dfs
    length_alignment_dfs = []
    for l in args.lengths:
        m_length = alignment_df['length'] == l
        l_df = pd.DataFrame(alignment_df[m_length])
        length_alignment_dfs.append((l, l_df))

    all_profiles_df = parallel.apply_parallel_iter(length_alignment_dfs, args.num_cpus, 
        get_metagene_profile, args, progress_bar=True)
    all_profiles_df = pd.concat(all_profiles_df)
    
    utils.write_df(all_profiles_df, args.out, index=False)

if __name__ == '__main__':
    main()

