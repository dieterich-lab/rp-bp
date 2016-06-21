#! /usr/bin/env python3

import argparse
import logging
import os
import sys
import time

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse

import pybedtools
import pysam
import tqdm

import misc.bio as bio
import misc.parallel as parallel
import misc.utils as utils

import misc.external_sparse_matrix_list as external_sparse_matrix_list

default_num_cpus = 2
default_num_orfs = 0
default_num_groups = 100
default_tmp = None # utils.abspath("tmp")

default_lengths = []
default_offsets = []

default_seqname_prefix = ''

def get_p_sites(bam_file, periodic_lengths, offsets):
    """ Given a bam file of mapped riboseq reads, this function filters
        out the reads of non-periodic length, adjusts the start and end
        positions based on strand, and then shifts the remaining reads
        based on the length-specific offset.
        
        Args:
            bam_file (string) : the path to the mapped riboseq reads
            
            periodic_lengths (list-like) : a list of lengths to keep
            
            offsets (list-like) : the distance to shift each read of the
                respective length. the order here must match that in
                periodic_lengths
                
        Returns:
            pd.DataFrame : a data frame containing the transformed reads,
                sorted by chrom and start
    """
    msg = "Reading BAM file"
    logging.info(msg)

    bam = pysam.AlignmentFile(bam_file)
    alignments = bam.fetch()
    num_alignments = bam.count()

    logging.info("Processing alignments")

    lengths = np.zeros(num_alignments, dtype=int)
    starts = np.zeros(num_alignments, dtype=int)
    ends = np.zeros(num_alignments, dtype=int)
    seqs = [""] * num_alignments
    strands = ["+"] * num_alignments

    for i, a in enumerate(tqdm.tqdm(alignments, leave=True, file=sys.stdout, total=num_alignments)):
        starts[i] = a.reference_start
        ends[i] = a.reference_end
        lengths[i] = a.qlen
        seqs[i] = a.reference_name

        if a.is_reverse:
            strands[i] = "-"

    # The data frame will later be converted to BED6, so put the fields in the
    # correct order.
    map_df = pd.DataFrame()
    map_df['seqname'] = seqs
    map_df['start'] = starts
    map_df['end'] = ends
    map_df['id'] = "."
    map_df['score'] = "."
    map_df['strand'] = strands
    map_df['length'] = lengths

    msg = "Filtering reads by length"
    logging.info(msg)
    
    # now, filter based on lengths
    m_length = map_df['length'].isin(periodic_lengths)
    map_df = map_df[m_length]

    # now, we need to update the starts and ends based on the strand
    msg = "Updating coordinates based on offsets"
    logging.info(msg)

    # if the strand is positive, the end is start+1
    # if the strand is negative, the start is end-1
    m_positive = map_df['strand'] == '+'
    m_negative = map_df['strand'] == '-'
    
    # first, shift in the appropriate direction
    for i in range(len(periodic_lengths)):
        length = periodic_lengths[i]
        offset = offsets[i]
        
        m_length = map_df['length'] == length
        
        # adjust the start of forward strand
        map_df.loc[m_positive & m_length, 'start'] = (
                map_df.loc[m_positive & m_length, 'start'] + offset)
        
        # adjust the ends of negative strand
        map_df.loc[m_negative & m_length, 'end'] = (
                map_df.loc[m_negative & m_length, 'end'] - offset)

    # finally, we only care about the 5' end of the read, so discard everything else
    msg = "Discarding 3' end of reads"
    logging.info(msg)
    
    map_df.loc[m_positive, 'end'] = map_df.loc[m_positive, 'start'] + 1
    map_df.loc[m_negative, 'start'] = map_df.loc[m_negative, 'end'] - 1

    # now sort everything
    msg = "Sorting reads by coordinates"
    logging.info(msg)
    
    map_df = map_df.sort_values(['seqname', 'start'])

    # and we only want the BED6 fields
    map_df = map_df[bio.bed6_field_names]
    
    return map_df

def get_orf_profile(orf_group):
    """ Given a "coverage_group", which contains all of the relative positions of reads
        which intersect a particular ORF. The reads at each position are counted. The
        profile is flipped for negative-strand ORFs.

        Args:
            orf_group (pd.DataFrame) : all of the reads which map to one ORF.

        Returns:
            string : the id of this ORF
            scipy.sparse_matrix : a sparse representation of the profile for this ORF
    """
    orf_name = orf_group['orf_id'].iloc[0]

    # pull out the profile
    orf_len = np.max(orf_group['orf_len'])
    profile = np.zeros(orf_len+1)
    for p in orf_group['position']:
        profile[p] += 1

    # flip, if needed
    strand = orf_group['strand'].iloc[0]
    if strand == '-':
        profile = profile[::-1]

    # sparsify and return
    sparse_profile = scipy.sparse.csr_matrix(profile)
    return (orf_name, sparse_profile)

def get_exon_index_rel_pos(abs_pos, exon_starts, seg_starts):
    """ This function determines the relative position of a given genomic
        position, relative to the provided exon and segment starts.

        Args:
            abs_pos (int) : the genomic position
            exon_starts, seg_starts (numpy.arrays of ints) : the absolute
                (genomic) starts and relative starts, respectively

        Returns:
            int : the exon into which the genomic position falls
            int : the relative position
    """
    # first, find the correct segment (which is also the correct exon)
    exon_index = np.searchsorted(exon_starts, abs_pos, side="right") - 1
    
    # get the offset from the end of the segment
    delta = abs_pos - exon_starts[exon_index]
    
    # move from the "start" of the segment
    rel_pos = seg_starts[exon_index] + delta
    
    return (exon_index, rel_pos)

def get_aligned_relative_position(interval):
    """ Given a pybedtools.intersect "Interval" object, this function parses the
        genomic coordinates of the ORF (the "A" feature, interpreted as a BED12
        record). It then finds the relative position of the matching read (the "B"
        feature).

        Args:
            interval (pybedtools.Interval) : an intersection result. The "A" 
                feature is interpreted as a BED12 record.

        Returns:
            dictionary : with the fields: orf_id, position, strand
    """
    start = int(interval[1])
    exon_lengths = interval[10]
    exon_relative_starts = interval[11]

    exon_lengths = np.array([int(l) for l in exon_lengths.split(',')], dtype=int)
    exon_starts = np.array([int(s) for s in exon_relative_starts.split(',')], dtype=int)
    exon_starts += start

    seg_starts = np.zeros(len(exon_lengths), dtype=int)
    seg_starts[1:] = np.cumsum(exon_lengths)[:-1]

    read_start = int(interval[13])
    
    exon_index, rel_pos = get_exon_index_rel_pos(read_start, exon_starts, seg_starts)
    
    orf_id = interval[3]
    strand = interval[5]
    
    ret = {
        'orf_id': orf_id,
        'position': rel_pos,
        'strand': strand
    }
    
    return ret

def get_orf_profiles(orfs, p_sites_file):
    start = time.perf_counter()
    
    # first, convert the ORFs into BED12
    orfs_bed = pybedtools.BedTool.from_dataframe(orfs[bio.bed12_field_names])
    
    # lazily get the intersection
    intersect_stream = orfs_bed.intersect(p_sites_file, split=True, s=True, stream=True, wo=True)
    aligned_positions = list(map(get_aligned_relative_position, intersect_stream))
    creation_time = time.perf_counter() - start

    msg = "Time creating relative alignments: {}".format(creation_time)
    logging.debug(msg)

    #pybedtools.helpers.close_or_delete(split_orfs)

    logging.debug("orfs_bed: {}".format(orfs_bed.fn))

    start = time.perf_counter()
    msg = "Number of alignments: {}".format(len(aligned_positions))
    logging.debug(msg)
    if len(aligned_positions) == 0:
        return None
    aligned_positions_df = pd.DataFrame(aligned_positions)
    
    convert_time = time.perf_counter() - start
    msg = "Time in conversion: {}".format(convert_time)
    logging.debug(msg)

    # we need to add the lengths
    aligned_positions_df = aligned_positions_df.merge(orfs[ ['id', 'orf_len'] ], left_on='orf_id', right_on='id')

    # pull out the profiles
    orf_groups = aligned_positions_df.groupby('orf_id')

    start = time.perf_counter()
    profiles = orf_groups.apply(get_orf_profile)
    profile_time = time.perf_counter() - start
    msg = "Time in profile call: {}".format(profile_time)
    logging.debug(msg)

    return profiles

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script uses bedtools to construct the profile for each ORF. It "
        "first adjusts the mapped read positions to properly align with the P-sites. Second, "
        "it uses bedtools to find the coverage of each position in each exon of each ORF. "
        "Finally, the ORF exons are glued together to find the profile of the entire ORF.")
    
    parser.add_argument('bam', help="The bam file including filtered (unique, etc.) alignments")
    parser.add_argument('orfs', help="The (bed12) file containing the ORFs")
    parser.add_argument('out', help="The (mtx) output file containing the ORF profiles")

    parser.add_argument('-l', '--lengths', help="If any values are given, then only reads "
        "which have those lengths will be included in the signal construction.", type=int,
        default=default_lengths, nargs='*')
    parser.add_argument('-o', '--offsets', help="The 5' end of reads will be shifted by this "
        "amount. There must be one offset value for each length (given by the --lengths "
        "argument.", type=int, default=default_offsets, nargs='*')
       
    parser.add_argument('-p', '--num-cpus', help="The number of processes to use", 
        type=int, default=default_num_cpus)
    parser.add_argument('-k', '--num-orfs', help="If  n>0, then only the first n orfs "
        "will be processed.", type=int, default=default_num_orfs)
    parser.add_argument('-g', '--num-groups', help="The number of groups into which to split "
        "the ORFs. More groups means the progress bar is updated more frequently but incurs "
        "more overhead because of the parallel calls.", type=int, default=default_num_groups)

    parser.add_argument('--seqname-prefix', help="If present, this string will be prepended "
        "to the seqname field of the ORFs.", default=default_seqname_prefix)

    parser.add_argument('--tmp', help="The temp directory for pybedtools", default=default_tmp)
        
    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    programs = ['bamToBed']
    utils.check_programs_exist(programs)

    if args.tmp is not None:
        os.makedirs(args.tmp, exist_ok=True)
        pybedtools.helpers.set_tempdir(args.tmp)
    
    # make sure the number of lengths and offsets match
    if len(args.lengths) != len(args.offsets):
        msg = "The number of --lengths and --offsets do not match."
        raise ValueError(msg)

    p_sites = get_p_sites(args.bam, args.lengths, args.offsets)

    msg = "Converting P-sites data frame to bed"
    logging.info(msg)

    p_sites_bed = pybedtools.BedTool.from_dataframe(p_sites)
    logging.debug("p_sites: {}".format(p_sites_bed.fn))
    #input("Go look at files!")

    # we do not need the data frame anymore, so save some memory

    msg = "Reading ORFs"
    logging.info(msg)

    orfs = bio.read_bed(args.orfs)
    orfs['seqname'] = orfs['seqname'].astype(str)

    if len(args.seqname_prefix) > 0:
        orfs['seqname'] = args.seqname_prefix + orfs['seqname']

    if args.num_orfs > 0:
        m_num_orfs = orfs['orf_num'] < args.num_orfs
        orfs = orfs[m_num_orfs]
     
    max_orf_len = np.max(orfs['orf_len'])
    msg = "The length of the longest ORF is: {}".format(max_orf_len)
    logging.debug(msg)

    # use bedtools to get the profiles
    msg = "Extracting profiles"
    logging.info(msg)

    profiles = parallel.apply_parallel_split(orfs, args.num_cpus, get_orf_profiles,
        p_sites_bed.fn, num_groups=args.num_groups, progress_bar=True)

    pybedtools.helpers.cleanup(verbose=True, remove_all=True)

    profiles = [p for p in profiles if p is not None]
    profiles = utils.flatten_lists(profiles)

    # get the largest transcript num
    logging.info("Copying profiles to matrix list")
    max_transcript_num = orfs['orf_num'].max()
    all_profiles = external_sparse_matrix_list.ExternalSparseMatrixList(max_transcript_num + 1)

    name_num_map = {name: int(num) for name,num in zip(orfs['id'], orfs['orf_num'])}

    for p in profiles:
        index = name_num_map[p[0]]
        sparse_profile = p[1]
        all_profiles[index] = sparse_profile

    # convert back to a single sparse matrix
    logging.info("Converting to single sparse matrix")

    # make sure we have enough columns to accomodate the longest ORF, even if
    # we do not have any reads from its batch.
    min_cols = max_orf_len + 1
    all_profiles = all_profiles.to_sparse_matrix(min_cols=min_cols)

    logging.info("Writing sparse matrix to disk")
    scipy.io.mmwrite(args.out, all_profiles)

if __name__ == '__main__':
    main()
