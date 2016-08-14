#! /usr/bin/env python3

import argparse
import collections
import itertools
import logging
import re
import numpy as np
import pandas as pd

import pyfasta
import pybedtools

import misc.bio as bio
import misc.bio_utils.gffread_utils as gffread_utils
import misc.logging_utils as logging_utils
import misc.parallel as parallel
import misc.utils as utils

logger = logging.getLogger(__name__)

default_novel_id_re = ""

default_start_codons = ['ATG']
default_stop_codons = ['TAA', 'TGA', 'TAG']

default_num_cpus = 1
default_num_transcripts = 0

def get_reverse_exon_index_gen_pos(rel_pos, exon_starts, exon_ends, seg_starts, seg_ends):
    """ This function converts a position in relative (i.e., transcript) coordinate into
        an absolute (i.e., genomic) coordinate. Furthermore, it identifies the exon to
        which the position belongs.

        !! This function assumes the exons and segments come from a transcript on the
            reverse strand !!

        Args:
            rel_pos (int) : the coordinate, in relative (transcript) space

            exon_starts, exon_ends (np.arrays of ints) : the absolute (genomic)
                coordinates of the exons

            seg_starts, seg_ends (np.arrays of ints) : the relative (transcript)
                coordinates of the exons

        Returns:
            int : the index of the exon into which rel_pos falls

            int : the abslute (genomic) position of rel_pos

        Example usage:
            exon, gen_pos = 
                get_reverse_exon_index_gen_pos(rel_start, exon_starts, exon_ends, seg_starts, seg_ends)
            
    """
    # first, find the correct segment
    seg_index = np.searchsorted(seg_starts, rel_pos, side="right") - 1

    # get the offset from the end of the segment
    delta_prime = seg_ends[seg_index] - rel_pos
    
    # flip to get the correct exon index
    exon_index = len(exon_starts) - seg_index - 1
    
    # move from the "start" of that exon
    gen_pos = exon_starts[exon_index] + delta_prime
    
    return (exon_index, gen_pos)

def get_reverse_orf_info(rel_start, rel_end, header):
    """ This function converts a relative start and end (i.e., relative ORF coordiantes)
        into a list of absolute coordinates. Additionally, it extracts the length and 
        relative genomic offset of each of the exons of the ORFs. Presumably, these are 
        used to create a BED12 file.

        !! This function assumes the ORF is on the reverse strand !!

        Args:
            rel_start (int) : the relative starting coordinate of the ORF

            rel_end (int) : the relative end coordinate of the ORF

            header (fasta_header) : a namedtuple containing at least the following fields:

                transcript_id, seqname, strand (strings) : the relevant information about
                    the transcript from which this ORF was extracted

                exon_starts, exon_ends (np.arrays of ints) : the absolute (genomic)
                    coordinates of the exons

                seg_starts, seg_ends (np.arrays of ints) : the relative (transcript)
                    coordinates of the exons

        Returns:
            pd.Series : a pandas series containing the ORF as a BED12 record            
            
    """
    start_exon, start_gen_pos = (
        get_reverse_exon_index_gen_pos(rel_start, header.exon_starts, header.exon_ends, 
            header.seg_starts, header.seg_ends))
    
    end_exon, end_gen_pos = (
        get_reverse_exon_index_gen_pos(rel_end, header.exon_starts, header.exon_ends, 
            header.seg_starts, header.seg_ends))
    
    num_exons = start_exon - end_exon + 1
    gen_starts = np.zeros(num_exons, dtype=int)
    gen_ends = np.zeros(num_exons, dtype=int)

    gen_starts[0] = end_gen_pos
    gen_ends[-1] = start_gen_pos
    
    if num_exons > 1:
        gen_ends[0] = header.exon_ends[end_exon]
        gen_starts[-1] = header.exon_starts[start_exon]
    
    for exon_index in range(end_exon+1, start_exon):
        s = header.exon_starts[exon_index]
        e = header.exon_ends[exon_index]

        gen_starts[exon_index-end_exon] = s
        gen_ends[exon_index-end_exon] = e
        
    # now, account for the base-0, half-open approach in BED
    gen_starts = gen_starts - 1
    
    # and collect the values we actually need for bed

    # we do not add 1 because of the half-open BED intervals
    exon_lengths = gen_ends - gen_starts
    exon_lengths = ','.join(str(e) for e in exon_lengths)
    
    rel_gen_starts = gen_starts - gen_starts[0]
    rel_gen_starts = ','.join(str(s) for s in rel_gen_starts)
    
    # use Mackowiak-stype orf_ids
    orf_id = "{}_{}:{}-{}:{}".format(header.transcript_id, header.seqname,
        gen_starts[0], gen_ends[-1], header.strand)
    
    ret = {
        "seqname": header.seqname,
        "start": gen_starts[0],
        "end": gen_ends[-1],
        "id": orf_id,
        "score": "0",
        "strand": header.strand,
        "thick_start": gen_starts[0],
        "thick_end": gen_ends[-1],
        "color": "0",
        "num_exons": num_exons,
        "exon_lengths": exon_lengths,
        "exon_genomic_relative_starts": rel_gen_starts,
        "cds_start": header.cds_start,
        "cds_end": header.cds_end
    }
    
    return pd.Series(ret)

def get_forward_exon_index_gen_pos(rel_pos, exon_starts, exon_ends, seg_starts, seg_ends):
    """ This function converts a position in relative (i.e., transcript) coordinate into
        an absolute (i.e., genomic) coordinate. Furthermore, it identifies the exon to
        which the position belongs.

        !! This function assumes the exons and segments come from a transcript on the
            forward strand !!

        Args:
            rel_pos (int) : the coordinate, in relative (transcript) space

            exon_starts, exon_ends (np.arrays of ints) : the absolute (genomic)
                coordinates of the exons

            seg_starts, seg_ends (np.arrays of ints) : the relative (transcript)
                coordinates of the exons

        Returns:
            int : the index of the exon into which rel_pos falls

            int : the abslute (genomic) position of rel_pos

        Example usage:
            exon, gen_pos = 
                get_reverse_exon_index_gen_pos(rel_start, exon_starts, exon_ends, seg_starts, seg_ends)
            
    """

    # first, find the correct segment (which is also the correct exon)
    seg_index = np.searchsorted(seg_starts, rel_pos, side="right") - 1    
    exon_index = seg_index
    
    # get the offset from the end of the segment
    delta = rel_pos - seg_starts[seg_index]
    
    # move from the "start" of the exon
    gen_pos = exon_starts[exon_index] + delta
    
    return (exon_index, gen_pos)

def get_forward_orf_info(rel_start, rel_end, header):
    """ This function converts a relative start and end (i.e., relative ORF coordiantes)
        into a list of absolute coordinates. Additionally, it extracts the length and 
        relative genomic offset of each of the exons of the ORFs. Presumably, these are 
        used to create a BED12 file.

        !! This function assumes the ORF is on the forward strand !!

        Args:
            rel_start (int) : the relative starting coordinate of the ORF

            rel_end (int) : the relative end coordinate of the ORF

            header (fasta_header) : a namedtuple containing at least the following fields:

                transcript_id, seqname, strand (strings) : the relevant information about
                    the transcript from which this ORF was extracted

                exon_starts, exon_ends (np.arrays of ints) : the absolute (genomic)
                    coordinates of the exons

                seg_starts, seg_ends (np.arrays of ints) : the relative (transcript)
                    coordinates of the exons

        Returns:
            pd.Series : a pandas series containing the ORF as a BED12 record            
    """

    start_exon, start_gen_pos = (
        get_forward_exon_index_gen_pos(rel_start, header.exon_starts, header.exon_ends, 
            header.seg_starts, header.seg_ends))
    
    end_exon, end_gen_pos = (
        get_forward_exon_index_gen_pos(rel_end, header.exon_starts, header.exon_ends, 
            header.seg_starts, header.seg_ends))
    
    num_exons = end_exon - start_exon + 1
    gen_starts = np.zeros(num_exons, dtype=int)
    gen_ends = np.zeros(num_exons, dtype=int)

    # first, take care of the end points
    gen_starts[0] = start_gen_pos
    gen_ends[-1] = end_gen_pos
    
    # now, if we had more than one exon, take care of the other sides of the first and last exon
    if num_exons > 1:
        gen_ends[0] = header.exon_ends[start_exon]
        gen_starts[-1] = header.exon_starts[end_exon]
    
    # and just copy over the exons in the middle
    for exon_index in range(start_exon+1, end_exon):
        s = header.exon_starts[exon_index]
        e = header.exon_ends[exon_index]

        gen_starts[exon_index-start_exon] = s
        gen_ends[exon_index-start_exon] = e

    # now, account for the base-0, half-open approach in BED
    gen_starts = gen_starts - 1
        
    # and collect the values we actually need for bed

    # we do not add 1 because of the half-open BED intervals
    exon_lengths = gen_ends - gen_starts
    exon_lengths = ','.join(str(e) for e in exon_lengths)
    
    rel_gen_starts = gen_starts - gen_starts[0]
    rel_gen_starts = ','.join(str(s) for s in rel_gen_starts)
    
    # use Mackowiak-stype orf_ids
    orf_id = "{}_{}:{}-{}:{}".format(header.transcript_id, header.seqname,
        gen_starts[0], gen_ends[-1], header.strand)
    
    ret = {
        "seqname": header.seqname,
        "start": gen_starts[0],
        "end": gen_ends[-1],
        "id": orf_id,
        "score": ".",
        "strand": header.strand,
        "thick_start": gen_starts[0],
        "thick_end": gen_ends[-1],
        "color": ".",
        "num_exons": num_exons,
        "exon_lengths": exon_lengths,
        "exon_genomic_relative_starts": rel_gen_starts,
        "cds_start": header.cds_start,
        "cds_end": header.cds_end
    }
    
    return pd.Series(ret)

def get_orf_end(start, stop_pos):
    """ This function finds the position of the first downstream, in-frame 
        stop for the given start codon. It then corrects this position because 
        BED files _do not_ include any of the stop codon in the annotation, and
        the coordinate ranges are treated as base-0, half-open. For example:
    
        a start position of 10 and end position of 20 covers 10 bases, [10,20),

        this interval starts at the *eleventh* base in the sequence

        Args:
            start (int) : the relative start position of the ORF. That is, the
                relative position of the "N" in "NUG"

            stop_pos (np.array of ints) : the relative positions of all stop
                codons. For example, the relative position of the "T" in "TAA".

        Returns:
            int: the relative position of the last base before the first
                downstream, in-frame stop codon

            OR None, if there are no downstream, in-frame stops
                
        
    """
    diff = stop_pos - start
    is_inframe_stop = (diff > 0) & (diff % 3 == 0)
    matches = np.where(is_inframe_stop==True)[0]
    
    if len(matches) == 0:
        return None
    
    # We need to subtract 1 here. I am still working out why, but it works.
    return (stop_pos[matches[0]] - 1)

def get_orfs_rel_pos(seq, start_codons_re, stop_codons_re):
    """ This function extracts the relative position of all ORFs from the given
        sequence. It assumes the sequence does not include any whitespace, and
        that the regular expressions properly identify start and stop codons.
        For example, if seq has already been transcribed, then the start codon
        should be "AUG".

        Args:
            seq (string) : the (untranslated) sequence

            start_codons_re, stop_codons_re (compiled regular expression):
                regular expressions which identify start and stop codons, respectively

        Returns:
            2d np.array of ints: the rows correspond to the 

        Example usage:
            orf_pos = get_orfs_rel_pos(seq, start_codons_re, stop_codons_re)

            first_orf_start = orf_starts[0, 0]
            first_orf_end = orf_ends[0, 1]
    """
    # add one so we consistently work in base-1
    
    # BED is in base-0, half-open, so add 1 to stop but not to start
    # however, our coordinates ARE in base-1, so work with those for the time being
    start_pos = np.array([m.start() for m in start_codons_re.finditer(seq)]) + 1
    stop_pos = np.array([m.start() for m in stop_codons_re.finditer(seq)]) + 1
    
    # pull out the matching ends for each start
    orf_ends = [get_orf_end(s, stop_pos) for s in start_pos]

    # zip them together
    orfs = zip(start_pos, orf_ends)

    # remove the starts which do not have a downstream stop
    filtered_orfs = list(filter(lambda x: x[1] is not None, orfs))
    return np.asarray(filtered_orfs)

def extract_orfs(header_seq, start_codons_re, stop_codons_re, ignore_parsing_errors):
    header, seq = header_seq
    seq = str(seq)

    try:
        h = gffread_utils.parse_header(header)
    except ValueError as ve:
        if ignore_parsing_errors:
            msg = "Could not parse header: '{}'".format(header)
            logger.warning(msg)
            return None
        else:
            raise
        
    orf_pos = get_orfs_rel_pos(seq, start_codons_re, stop_codons_re)

    # use function pointers to pull out the correct orf information based on strand
    if h.strand == '+':
        get_orf_info = get_forward_orf_info
    else:
        get_orf_info = get_reverse_orf_info

    orf_info = [get_orf_info(o[0], o[1], h) 
                for o in orf_pos]

    orf_info = pd.DataFrame(orf_info)
    return orf_info



def extract_canonical_orf(header_seq, start_codons_re, stop_codons_re, ignore_parsing_errors):
    """ This function checks if the given header includes an annotated CDS. If so,
        it extracts the information for that ORF.

        Args:
            header_seq (tuple) : the first item of the tuple is the fasta header
                (without the ">"), and the second item is the pyfasta sequence object.
                The argument is given this way to make it easy to use this function
                while iterating through the pyfasta file.

            start_codons_re, stop_codons_re (compiled re) : compiled reg_ex
                objects used to locate start and stop codons, respectively

        Returns:
            pd.Series : a pandas series containing the ORF as a BED12 record
                
                OR

                None, if the header does not specific a canonical ORF
    """
    header, seq = header_seq
    seq = str(seq)

    
    try:
        h = gffread_utils.parse_header(header)
    except ValueError as ve:
        if ignore_parsing_errors:
            msg = "Could not parse header: '{}'".format(header)
            logger.warning(msg)
            return None
        else:
            raise

    if h.cds_start is None:
        return None

    # we will still need the stop codon positions, so grab them
    orf_pos = get_orfs_rel_pos(seq, start_codons_re, stop_codons_re)

    # now, find the ORF which corresponds to the canonical CDS
    canonical_orf = np.where(orf_pos == h.cds_start)
    if len(canonical_orf[0]) == 0:
        possible_start = seq[h.cds_start: h.cds_start+3]
        msg = ("The canoncial ORF could not be found. This could be caused, for example, "
                "if the canonical ORF uses an alternative start codon. The codon at the "
                "annotated start position is '{}'. Header: '{}'".format(possible_start, 
                header))
        logger.error(msg)
        return None

    # and pull the actual index out of the np.where result
    matching_rows = canonical_orf[0]
    matching_columns = canonical_orf[1]
    
    # pick out the row where the !start! matches the cds_start
    r = np.where(matching_columns == 0)[0]
    
    if len(r) == 0:
        possible_start = seq[h.cds_start: h.cds_start+3]
        msg = ("The canonical ORF could not be found. This could be caused, for example, "
                "if the canonical ORF uses an alternative start codon. The codon at the "
                "annotated start position is '{}'. Header: '{}'".format(possible_start, 
                header))
        logger.error(msg)
        return None

    r = r[0]
    canonical_orf_index = matching_rows[r]
    canonical_orf = orf_pos[canonical_orf_index]

    # use function pointers to pull out the correct orf information based on strand
    if h.strand == '+':
        get_orf_info = get_forward_orf_info
    else:
        get_orf_info = get_reverse_orf_info

    canonical_orf_info = get_orf_info(canonical_orf[0], canonical_orf[1], h)
    return canonical_orf_info


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script extracts all of the ORFs from the provided transcript "
        "fasta file into a BED file. Furthermore, it labels each ORF with respect to the "
        "annotated canonical coding sequences.")
    parser.add_argument('transcript_fasta', help="The transcript fasta file produced by "
        "gffread (from the Cufflinks package)")
    parser.add_argument('out', help="The output (bed.gz) file")

    parser.add_argument('--start-codons', help="The sequences treated as start codons",
        nargs='+', default=default_start_codons)
    parser.add_argument('--stop-codons', help="The sequences treated as stop codons",
        nargs='+', default=default_stop_codons)

    parser.add_argument('--ignore-parsing-errors', help="If this flag is present, then "
        "headers which do not parse will be skipped. A warning is printed for each "
        "header skipped in this manner.", action='store_true')

    parser.add_argument('--novel-id-re', help="This option is a regular expression for "
        "the identifiers of novel ORFs. ORFs are annotated as novel if their id's match "
        "this expression and they would otherwise be annotated as 'noncoding' or "
        "'suspect_overlap'. If no expression is given, no ORFs are annotated as novel.", 
        default=default_novel_id_re)

    parser.add_argument('--num-cpus', help="The number of CPUs to use", type=int, 
        default=default_num_cpus)

    parser.add_argument('--num-transcripts', help="If n>0, then only the first n "
        "transcripts will be processed.", type=int, default=default_num_transcripts)

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    programs = ['intersectBed']
    utils.check_programs_exist(programs)

    msg = "Compiling the start and stop codon regular expressions"
    logger.info(msg)
    
    start_codons_re = '|'.join(args.start_codons)
    stop_codons_re = '|'.join(args.stop_codons)

    start_codons_re = re.compile(start_codons_re)
    stop_codons_re = re.compile(stop_codons_re)

    msg = "Reading the transcript sequences"
    logger.info(msg)

    #transcript_fasta = pyfasta.Fasta(args.transcript_fasta)
    transcript_fasta = bio.get_read_iterator(args.transcript_fasta)

    # check if we only want to process a few of the transcripts

    # we must use a list; otherwise, we will not have the same order
    # for transcripts later
    if args.num_transcripts > 0:
        transcripts = list(itertools.islice(transcript_fasta, args.num_transcripts))
        total = args.num_transcripts
    else:
        transcripts = list(transcript_fasta)
        total = len(transcripts)

    msg = "Extracting ORFs"
    logger.info(msg)

    orfs = parallel.apply_parallel_iter(transcripts, args.num_cpus, extract_orfs,
                                     start_codons_re, stop_codons_re,
                                     args.ignore_parsing_errors,
                                     progress_bar=True, total=total)

    # filter out any Nones
    # these could come up if we ignore parsing errors
    orfs = [o for o in orfs if o is not None]

    msg = "Combining ORFs into large data frame"
    logger.info(msg)
    orfs = pd.concat(orfs)

    msg = "Sorting ORFs"
    logger.info(msg)
    orfs = orfs.sort_values(['seqname', 'start'])

    msg = "Removing duplicate ORFs"
    logger.info(msg)
    duplicate_fields = ['seqname', 'start', 'end', 'strand', 'num_exons']
    orfs = orfs.drop_duplicates(subset=duplicate_fields)

    msg = "Found {} unique ORFs".format(len(orfs))
    logger.info(msg)

    msg = "Extracting canonical ORFs"
    logger.info(msg)

    # extract the canonical ORFs based on the annotation
    canonical_orfs = parallel.apply_parallel_iter(transcripts, args.num_cpus,
                                               extract_canonical_orf, 
                                               start_codons_re, stop_codons_re,
                                               args.ignore_parsing_errors,
                                               progress_bar=True, total=total)

    canonical_orfs = [o for o in canonical_orfs if o is not None]
    canonical_orfs = pd.DataFrame(canonical_orfs)
    canonical_orfs = canonical_orfs.sort_values(['seqname', 'start'])

    # update the field names so we can distinguish between the ORFs after bedutils/join
    canonical_orfs = (canonical_orfs[bio.bed12_field_names]).copy()
    canonical_orfs.columns = ['canonical_{}'.format(c) for c in bio.bed12_field_names]
    canonical_orfs_bed = pybedtools.BedTool.from_dataframe(canonical_orfs)

    # and extract field names for parsing
    intersect_fields = bio.bed12_field_names + list(canonical_orfs.columns) + ['overlap']

    orfs['orf_type'] = None

    
    ### canonical

    msg = "Labeling canonical, annotated CDSs"
    logger.info(msg)

    # pull out the unannotated orfs
    m_remaining = orfs['orf_type'].isnull()
    orfs_bed = pybedtools.BedTool.from_dataframe(orfs.loc[m_remaining, bio.bed12_field_names])

    msg = "Number of BED records: {}".format(len(orfs_bed))
    logger.debug(msg)

    # intersect with the canonical ORFs
    # split: include BED12 information (i.e., split the exons)
    # s: require matches be on the same strand
    # wo: write both entries, as well as their overlap
    # f 1.0 : require coverage of the entire "A" feature
    # F 1.0 : require coverage of the entire "B" feature
    i_canonical = orfs_bed.intersect(canonical_orfs_bed, split=True, s=True, wo=True, f=1.0, F=1.0)
    i_canonical_df = i_canonical.to_dataframe(names=intersect_fields)

    pybedtools.helpers.close_or_delete(i_canonical)

    msg = "Found {} intersecting records".format(len(i_canonical_df))
    logger.debug(msg)
    
    # pull out the orf_ids of everything with a match
    i_ids = i_canonical_df['id']
    m_canonical = orfs['id'].isin(i_ids)

    # and label those orfs
    orfs.loc[m_canonical, 'orf_type'] = 'canonical'

    msg = "Found {} canonical ORFs".format(sum(m_canonical))
    logger.info(msg)
    
    ### canonical_extended
    
    msg = "Finding \"extended\" canonical ORFs, which completely cover annotated CDSs"
    logger.info(msg)

    # pull out the remaining unannotated orfs
    m_remaining = orfs['orf_type'].isnull()
    orfs_bed = pybedtools.BedTool.from_dataframe(orfs.loc[m_remaining, bio.bed12_field_names])
    
    msg = "Number of BED records: {}".format(len(orfs_bed))
    logger.debug(msg)

    # intersect with canonical ORFs
    # split: include BED12 information (i.e., split the exons)
    # s: require matches be on the same strand
    # wo: write both entries, as well as their overlap
    # F 1.0 : require coverage of the entire "B" feature
    # 
    # Since we have already pulled out the ORFs which _exactly_ map to canonical CDSs,
    # the remaining ORFs which completely cover a canonical orf must _extend_ it in
    # some direction.
    #
    # In theory, this could be a huge, out-of-frame ORF relative to the annotation;
    # however, we do not isolate that case here.
    i_canonical_extended = orfs_bed.intersect(canonical_orfs_bed, split=True, s=True, wo=True, F=1.0)
    i_canonical_extended_df = i_canonical_extended.to_dataframe(names=intersect_fields)

    pybedtools.helpers.close_or_delete(i_canonical_extended)

    msg = "Found {} intersecting records".format(len(i_canonical_extended_df))
    logger.debug(msg)

    # pull out the orf_ids of everything with a match
    i_ids = i_canonical_extended_df['id']
    m_canonical_extended = orfs['id'].isin(i_ids)

    # and label them
    orfs.loc[m_canonical_extended, 'orf_type'] = 'canonical_extended'

    msg = "Found {} canonical_extended ORFs".format(sum(m_canonical_extended))
    logger.info(msg)

    ### canonical_truncated
    msg = "Finding \"canonical_truncated\" and \"within\" ORFs"
    logger.info(msg)

    # pull out the remaining unannotated orfs
    m_remaining = orfs['orf_type'].isnull()
    orfs_bed = pybedtools.BedTool.from_dataframe(orfs.loc[m_remaining, bio.bed12_field_names])

    msg = "Number of BED records: {}".format(len(orfs_bed))
    logger.debug(msg)

    # intersect with the canonical ORFs
    # split: include BED12 information (i.e., split the exons)
    # s: require matches be on the same strand
    # wo: write both entries, as well as their overlap
    # f 1.0 : require coverage of the entire "A" feature
    #
    # Since we have already pulled out the ORFs which _exactly_ map to canonical CDSs,
    # the remaining ORFs which are completely covered by a canonical CDS are either
    # truncated versions of that transcript or out-of-frame ORFs within it.
    i_canonical_truncated = orfs_bed.intersect(canonical_orfs_bed, split=True, s=True, wo=True, f=1.0)
    i_canonical_truncated_df = i_canonical_truncated.to_dataframe(names=intersect_fields)

    pybedtools.helpers.close_or_delete(i_canonical_truncated)

    msg = "Found {} intersecting records".format(len(i_canonical_truncated_df))
    logger.debug(msg)

    # truncated ORFs on the positive strand have the same end as the canonical ORF
    # while those on the negative strand have the same "start" (really end, but 
    # "start" because of BED notation)
    m_matching_start = i_canonical_truncated_df['start'] == i_canonical_truncated_df['canonical_start']
    m_matching_end = i_canonical_truncated_df['end'] == i_canonical_truncated_df['canonical_end']

    m_forward = i_canonical_truncated_df['strand'] == '+'
    m_reverse = i_canonical_truncated_df['strand'] == '-'
    m_truncated_forward = m_forward & m_matching_end

    m_truncated_reverse = m_reverse & m_matching_start
    m_canonical_truncated = m_truncated_forward | m_truncated_reverse

    # everything left is completely overlapped by some canonical ORF, but does not share a start or end

    # label "within" first. that way, if an ORF is associated with multiple canonical CDS regions, then
    # "truncated" takes precedence
    m_within = ~m_canonical_truncated
    i_ids = i_canonical_truncated_df.loc[m_within, 'id']
    m_within = orfs['id'].isin(i_ids)
    orfs.loc[m_within, 'orf_type'] = 'within'

    i_ids = i_canonical_truncated_df.loc[m_canonical_truncated, 'id']
    m_ct = orfs['id'].isin(i_ids)
    orfs.loc[m_ct, 'orf_type'] = 'canonical_truncated'

    msg = "Found {} \"within\" ORFs".format(sum(m_within & ~m_canonical_truncated))
    logger.info(msg)

    msg = "Found {} \"canonical_truncated\" ORFs".format(sum(m_canonical_truncated))
    logger.info(msg)

    ### leader/trailer overlaps
    msg = "Finding out-of-frame ORFs which overlap canonical CDSs"
    logger.info(msg)

    # pull out the remaining orfs
    m_remaining = orfs['orf_type'].isnull()
    orfs_bed = pybedtools.BedTool.from_dataframe(orfs.loc[m_remaining, bio.bed12_field_names])
    
    msg = "Number of BED records: {}".format(len(orfs_bed))
    logger.debug(msg)

    # find ORFs which have some overlap with the canonical CDS regions
    i_utr_overlap = orfs_bed.intersect(canonical_orfs_bed, split=True, s=True, wo=True)
    i_utr_overlap_df = i_utr_overlap.to_dataframe(names=intersect_fields)

    pybedtools.helpers.close_or_delete(i_utr_overlap)

    msg = "Found {} intersecting records".format(len(i_utr_overlap_df))
    logger.debug(msg)

    # here, we have 4 cases:
    #   forward strand, 5' overlap : start < canonical_start
    #   forward strand, 3' overlap : end > canonical_end
    #   reverse strand, 5' overlap : end > canonical_end
    #   reverse strand, 3' overlap : start < canonical_start
    m_forward = i_utr_overlap_df['strand'] == '+'
    m_reverse = i_utr_overlap_df['strand'] == '-'

    m_start = i_utr_overlap_df['start'] < i_utr_overlap_df['canonical_start']
    m_end = i_utr_overlap_df['end'] > i_utr_overlap_df['canonical_end']

    m_five_prime_overlap = (m_forward & m_start) | (m_reverse & m_end)
    m_three_prime_overlap = (m_forward & m_end) | (m_reverse & m_start)

    # and label the overlapping ORFs
    i_ids = i_utr_overlap_df.loc[m_five_prime_overlap, 'id']
    m_fp = orfs['id'].isin(i_ids)
    orfs.loc[m_fp, 'orf_type'] = 'five_prime_overlap'

    i_ids = i_utr_overlap_df.loc[m_three_prime_overlap, 'id']
    m_tp = orfs['id'].isin(i_ids)
    orfs.loc[m_tp, 'orf_type'] = 'three_prime_overlap'

    msg = "Found {} \"five_prime_overlap\" ORFs".format(sum(m_five_prime_overlap))
    logger.info(msg)

    msg = "Found {} \"three_prime_overlap\" ORFs".format(sum(m_three_prime_overlap))
    logger.info(msg)


    ### leader/trailer ORFs
    msg = "Finding ORFs entirely in the 5' leaders and 3' trailers, and noncoding ORFs"
    logger.info(msg)

    m_remaining = orfs['orf_type'].isnull()
    orfs_bed = pybedtools.BedTool.from_dataframe(orfs.loc[m_remaining, bio.bed12_field_names])

    msg = "Number of BED records: {}".format(len(orfs_bed))
    logger.debug(msg)

    # for closest to work, we have to split the "B" file (canonical orfs) into separate exons
    # I cannot find documentation for why...
    canonical_orfs_bed6 = canonical_orfs_bed.bed6()
    canonical_orfs_bed6_sorted = canonical_orfs_bed6.sort()
    canonical_bed6_fields = ['canonical_{}'.format(c) for c in bio.bed6_field_names]
    bed6_closest_fields = bio.bed12_field_names + canonical_bed6_fields + ['distance']

    # now, find the relationship of the remaining ORFs to the annotated CDS regions
    # based on those, label the ORFs as either 5' or 3'
    closest = orfs_bed.closest(canonical_orfs_bed6_sorted, s=True, D='b', t='first')
    closest_df = closest.to_dataframe(names=bed6_closest_fields, index_col=False)

    pybedtools.helpers.close_or_delete(closest)

    msg = "Found {} intersecting records".format(len(closest_df))
    logger.debug(msg)

    m_overlap = closest_df['distance'] == 0
    m_three_prime = closest_df['distance'] > 0
    m_five_prime = closest_df['distance'] < 0

    # additionally, check for noncoding ORFs which do not overlap coding ORFs
    m_no_cds = orfs['cds_start'].isnull()


    i_ids = closest_df.loc[m_five_prime, 'id']
    m_fp = orfs['id'].isin(i_ids)
    orfs.loc[m_fp & ~m_no_cds, 'orf_type'] = 'five_prime'

    i_ids = closest_df.loc[m_three_prime, 'id']
    m_tp = orfs['id'].isin(i_ids)
    orfs.loc[m_tp & ~m_no_cds, 'orf_type'] = 'three_prime'

    i_ids = closest_df.loc[m_overlap, 'id']
    m_o = orfs['id'].isin(i_ids)

    # if we specified a form for novel IDs, annotate them as such
    if len(args.novel_id_re) > 0:
        novel_id_re = re.compile(args.novel_id_re)
        m_novel = orfs['id'].str.match(novel_id_re)

        orfs.loc[m_novel & m_fp & m_no_cds, 'orf_type'] = 'novel'
        orfs.loc[m_novel & m_tp & m_no_cds, 'orf_type'] = 'novel'
        orfs.loc[m_novel & m_o, 'orf_type'] = 'novel_suspect_overlap'
    else:
        m_novel = False * len(orfs)
    
    orfs.loc[m_fp & m_no_cds & ~m_novel, 'orf_type'] = 'noncoding'
    orfs.loc[m_tp & m_no_cds & ~m_novel, 'orf_type'] = 'noncoding'
    orfs.loc[m_o & ~m_novel, 'orf_type'] = 'suspect_overlap'

    msg = "Found {} \"five_prime\" ORFs".format(sum(m_fp & ~m_no_cds))
    logger.info(msg)

    msg = "Found {} \"three_prime\" ORFs".format(sum(m_tp & ~m_no_cds))
    logger.info(msg)

    msg = "Found {} \"suspect_overlap\" ORFs".format(sum(m_o))
    logger.info(msg)

    num_noncoding = sum(m_tp & m_no_cds) + sum(m_fp & m_no_cds)
    msg = "Found {} \"noncoding\" ORFs".format(num_noncoding)
    logger.info(msg)


    m_remaining = orfs['orf_type'].isnull()
    msg = "{} ORFs are unlabeled".format(sum(m_remaining))
    logger.info(msg)

    msg = "Sorting labeled ORFs"
    logger.info(msg)

    orfs = orfs.sort_values(['seqname', 'start'])

    # now, get the lengths of all of the ORFs
    msg = "Getting ORF lengths"
    logger.info(msg)

    orf_lengths = parallel.apply_parallel(orfs, args.num_cpus, 
        bio.get_bed_12_feature_length, progress_bar=True)
    orf_lengths = np.array(orf_lengths) 
    orfs['orf_len'] = orf_lengths


    msg = "Writing ORFs as BED"
    logger.info(msg)

    # also, add a numeric index field
    orfs['orf_num'] = range(len(orfs))

    output_fields = bio.bed12_field_names + ['orf_type', 'orf_len', 'orf_num']
    bio.write_bed(orfs[output_fields], args.out)

if __name__ == '__main__':
    main()
