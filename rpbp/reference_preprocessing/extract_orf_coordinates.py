#! /usr/bin/env python3

import argparse
import collections
import numpy as np
import pandas as pd
import re

import bio_utils.bed_utils as bed_utils
import bio_utils.fastx_utils as fastx_utils
import misc.parallel as parallel
import misc.utils as utils

import misc.slurm as slurm
import logging
import misc.logging_utils as logging_utils
logger = logging.getLogger(__name__)

default_start_codons= ['ATG']
default_stop_codons = ['TAA', 'TGA', 'TAG']

# these fields will be used to detect duplicate ORFs
DUPLICATE_FIELDS = [
    'seqname', 
    'start', 
    'end', 
    'strand', 
    'num_exons', 
    'exon_lengths', 
    'exon_genomic_relative_starts'
]

orf_position = collections.namedtuple(
    "orf_position",
    "start,end"
)

def get_orf_positions(seq, start_codons_re, stop_codons_re):
    """ This function extracts the relative position of all ORFs from the given
        sequence. It assumes the sequence does not include any whitespace, and
        that the regular expressions properly identify start and stop codons.
        For example, if seq has already been transcribed, then the start codon
        should be "AUG".
        
        N.B. The ORFs *include* the first base in the start codon (e.g., "A" in
        "ATG") 

        Args:
            seq (string) : the (untranslated) sequence

            start_codons_re, stop_codons_re (compiled regular expression):
                regular expressions which identify start and stop codons, respectively

        Returns:
            list of orf_positions (a named 2-tuple with "start" and "end" fields)
            

        Example usage:
            orfs = get_orfs(seq, start_codons_re, stop_codons_re)

            first_orf_start = orf_starts[0, 0]
            first_orf_end = orf_ends[0, 1]
    """
    
    # these give the positions of the 
    start_pos = np.array([m.start() for m in start_codons_re.finditer(seq)])
    stop_pos = np.array([m.start() for m in stop_codons_re.finditer(seq)])
    
    # pull out the matching ends for each start
    orfs = [get_matching_stop_position(s, stop_pos) for s in start_pos]
    orfs = utils.remove_nones(orfs)
    return orfs

def get_matching_stop_position(start, stop_pos):
    """ This function finds the position of the first downstream, in-frame 
        stop for the given start codon. It returns the ORF indices as a tuple.
        
        N.B. The coordinates are given in bed-style *half-open* intervals, so
        the "start" base is *included* but the "stop" base is *excluded*.
        
        For example, "start" could point to the "A" in ATG while "stop" could
        point to the "T" in "TAA".

        Args:
            start (int) : the relative start position of the ORF. That is, the
                relative position of the "N" in "NUG"

            stop_pos (np.array of ints) : the relative positions of all stop
                codons. For example, the relative position of the "T" in "TAA".

        Returns:
            orf_position, which is a named 2-tuple with "start" and "end" fields.
                The positions are the (relative) position of the orf starting at start.

            OR None, if there are no downstream, in-frame stops
    """
    diff = stop_pos - start
    is_inframe_stop = (diff > 0) & (diff % 3 == 0)
    matches = np.where(is_inframe_stop==True)[0]
    
    if len(matches) == 0:
        return None
    
    return orf_position(start, stop_pos[matches[0]])

def get_orf_bed_entry(orf_position, transcript):
    """ This function takes the genomic start and end positions of an ORF, 
        and the transcript from which it originates, and returns a bed12+1
        entry corresponding just to the ORF.
        
        This function is largely a wrapper around bed_utils.retain_thick_only, 
        so please see that function for more details.
        
        Besided the normal bed12 fields, this function also includes 'orf_len',
        which gives the length of the ORF.
        
        N.B. This function *does not* take into strand information. That
        needs to be done before this function call.
        
        Args:
            orf_position (orf_position): the genomic start and end of the ORF
            
            transcript (dict-like):  a bed12+ entry which can be indexed as a
                dictionary (e.g., a dict or pd.Series). The following field 
                names must match those in bio.bed12_field_names:
                
                    start, end, thick_start, end, thick_end
                    num_exons, exon_lengths, exon_genomic_relative_starts
                    
                Additionally, the object must have a "copy" function.
                
        Returns:
            dict-like: a bed entry for the ORF
 
    """
    # make a copy
    orf = transcript.copy()
    
    # and trim everything except the ORF
    orf['thick_start'] = orf_position.start
    orf['thick_end'] = orf_position.end

    bed_utils.retain_thick_only(orf, inplace=True)
    
    orf_len = bed_utils.get_bed_12_feature_length(orf)
    orf['orf_len'] = orf_len
    
    # use Mackowiak-stype orf_ids
    orf_id = "{}_{}:{}-{}:{}".format(orf['id'], orf['seqname'],
        orf['start'], orf['end'], orf['strand'])
    orf['id'] = orf_id

    return orf

def get_orfs(transcript_and_sequence, start_codons_re, stop_codons_re):
    """ This function extracts all of the ORFs from the given transcript and
        sequence. It returns them as a bed12+1 data frame.
    """
    transcript, transcript_sequence = transcript_and_sequence
    transcript_length = len(transcript_sequence)
    
    # get the ORFs for this entry
    orf_rel_positions = get_orf_positions(
        transcript_sequence,
        start_codons_re,
        stop_codons_re
    )
    
    #if logger.isEnabledFor(logging.DEBUG):
    #   s = ["({},{})".format(o.start, o.end) for o in orf_rel_positions]
    #    s = ' '.join(s)
    #    logger.debug(s)

    # now, if the strand is negative, we need to "flip" the relative positions 
    if transcript['strand'] == '-':

        orf_rel_positions = [
            orf_position(
                start = transcript_length - o.end,
                end = transcript_length - o.start
            ) for o in orf_rel_positions
        ]
        

        # this is too much debug info; unit tests handle this
        #if logger.isEnabledFor(logging.DEBUG):
        #    s = ["({},{})".format(o.start, o.end) for o in orf_rel_positions]
        #    s = ' '.join(s)
        #    logger.debug(s)

    # we need the block information to convert between relative and genomic coordinates
    start = transcript['start']
    
    block_lengths = np.fromstring(
        transcript['exon_lengths'],
        sep=',',
        dtype=int
    )

    block_starts = np.zeros(len(block_lengths), dtype=int)
    block_starts[1:] = np.cumsum(block_lengths)[:-1]

    block_relative_starts = np.fromstring(
        transcript['exon_genomic_relative_starts'],
        sep=',',
        dtype=int
    )

    # for a discussion about why 
    # please see Issue #64: https://github.com/dieterich-lab/rp-bp/issues/64
    orf_gen_positions = [
        orf_position(
            start=bed_utils.get_gen_pos(
                o.start,
                start,
                block_lengths,
                block_starts,
                block_relative_starts),

            end=bed_utils.get_gen_pos(
                o.end-1,
                start,
                block_lengths,
                block_starts,
                block_relative_starts)+1
        ) for o in orf_rel_positions
    ]

    # construct a data frame from the record for each
    orfs = [get_orf_bed_entry(o, transcript) for o in orf_gen_positions]
    orfs = pd.DataFrame(orfs)
    
    return orfs

def get_transcript(transcript_id, transcripts_bed):
    """ This is a simple helper function to grab the right transcript out of
        the data frame for the iterator.
    """
    m_transcript = transcripts_bed['id'] == transcript_id
    transcript = transcripts_bed[m_transcript].iloc[0]
    return transcript

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script extracts all of the ORFs from the given transcripts. "
        "It writes the result as a bed12+1 file. The additional field, 'orf_len', gives "
        "the length of the respective ORF. It removes duplicate ORFs.\n\nN.B. The DEBUG "
        "output for this script is _very_ verbose. It is not recommended to run this "
        "script with that logging level.")

    parser.add_argument('transcripts_bed', help="The bed12 file containing the "
        "transcript information")

    parser.add_argument('transcripts_fasta', help="The fasta file containing the "
        "spliced transcript sequences")

    parser.add_argument('out', help="The output (bed12+1 gz) file")

    parser.add_argument('--start-codons', help="A list of codons which will be "
        "treated as start codons when extracting ORFs", nargs='+', 
        default=default_start_codons)

    parser.add_argument('--stop-codons', help="A list of codons which will be "
        "treated as stop codons when extracting ORFs", nargs='+', 
        default=default_stop_codons)

    slurm.add_sbatch_options(parser)
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    # check if we wanted to use slurm
    if args.use_slurm:
        cmd = ' '.join(sys.argv)
        slurm.check_sbatch(cmd, args=args)
        return


    msg = "Compiling start and stop codon regular expressions"
    logger.info(msg)
    
    start_codons_re = '|'.join(args.start_codons)
    stop_codons_re = '|'.join(args.stop_codons)

    start_codons_re = re.compile(start_codons_re)
    stop_codons_re = re.compile(stop_codons_re)

    msg = "Reading transcripts bed file"
    logger.info(msg)
    transcripts_bed = bed_utils.read_bed(args.transcripts_bed)
    
    msg = "Creating the sequence iterator"
    logger.info(msg)

    transcripts_fasta = fastx_utils.get_read_iterator(args.transcripts_fasta)

    transcripts_iter = ((get_transcript(transcript_header, transcripts_bed), transcript_sequence)
                    for (transcript_header, transcript_sequence) in transcripts_fasta)

    msg = "Finding all ORFs"
    logger.info(msg)

    orfs = parallel.apply_parallel_iter(transcripts_iter,
                                        args.num_cpus,
                                        get_orfs,
                                        start_codons_re, stop_codons_re,
                                        total = len(transcripts_bed),
                                        progress_bar=True)

    msg = "Joining ORFs in a large data frame"
    logger.info(msg)

    orfs = pd.concat(orfs)

    msg = "Removing duplicate ORFs"
    logger.info(msg)
    
    orfs = orfs.drop_duplicates(subset=DUPLICATE_FIELDS)

    msg = "Numbering remaining ORFs"
    logger.info(msg)

    orfs['orf_num'] = np.arange(len(orfs))

    msg = "Writing ORFs to disk"
    logger.info(msg)
    bed_utils.write_bed(orfs, args.out)

if __name__ == '__main__':
    main()
