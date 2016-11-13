#! /usr/bin/env python3

import argparse
import pandas as pd
import tqdm

import misc.bio_utils.bed_utils as bed_utils
import misc.logging_utils as logging_utils
import misc.parallel as parallel

import logging
import misc.logging_utils as logging_utils
logger = logging.getLogger(__name__)

default_num_cpus = 1

default_annotated_exons = None
default_label_prefix = ""
default_nonoverlapping_label = None

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script labels the ORFs found with extract-orf-coordinates "
        "based on their exon structure and relation to the annotated, canonical "
        "ORFs. It requires the exon blocks for the ORFs (created with "
        "split-bed12-blocks). It completely reads in the ORFs, so unless otherwise "
        "desired for some reason, the input and output files can be the same.")

    parser.add_argument('annotated_transcripts', help="The annotated transcripts "
        "for the genome, in bed12+ format")
    parser.add_argument('extracted_orfs', help="The ORFs extracted from the "
        "transcripts, in bed12+ format")
    parser.add_argument('orf_exons', help="The exon blocks for the ORFs, in "
        "bed6+ format")

    parser.add_argument('out', help="The output (bed12+.gz) file")

    parser.add_argument('-p', '--num-cpus', help="The number of CPUs to use for "
        "a few parts of the script", type=int, default=default_num_cpus)

    parser.add_argument('-f', '--filter', help="If this flag is given, then ORFs "
        "which are completely covered by an annotated transcript are discarded. "
        "Presumably, this is used to filter uninteresting ORFs from de novo "
        "assemblies.", action='store_true')

    parser.add_argument('-e', '--annotated-exons', help="If the --filter flag is "
        "given, the annotated transcript exons can optionally be provided with "
        "this option. If they are not given, they will be split from the annotated "
        "transcripts. That is generally not a very expensive operation relative to "
        "everything else in the labeling script. If --filter is not given, then "
        "these are ignored.", default=default_annotated_exons)

    parser.add_argument('-n', '--nonoverlapping-label', help="If this option is "
        "given, then ORFs which do not overlap the annotated transcripts at all "
        "will be given this label. Otherwise, they will be labeled as \"suspect\"",
        default=default_nonoverlapping_label)

    parser.add_argument('-l', '--label-prefix', help="This string is prepended "
        "to all labels assigned to ORFs. For example, it is a useful way to "
        "indicate ORFs from de novo assemblies are \"novel.\" In any case, this "
        "*is not* prepended to \"canonical\" ORFs.", default=default_label_prefix)
    
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "Reading annotated transcripts"
    logger.info(msg)
    annotated_transcripts = bed_utils.read_bed(args.annotated_transcripts)

    msg = "Reading extracted ORFs and exons"
    logger.info(msg)
    extracted_orfs = bed_utils.read_bed(args.extracted_orfs)
    extracted_orf_exons = bed_utils.read_bed(args.orf_exons)

    msg = "Found {} extracted ORFs with {} exons".format(len(extracted_orfs),
        len(extracted_orf_exons))
    logger.debug(msg)

    # check if we want to remove the extracted_orfs completely covered by
    # the annotated transcripts
    if args.filter:
        msg = ("Removing extracted ORFs which are completely covered by the "
            "annotated transcripts")
        logger.info(msg)
        
        # we need the annotated transcript exons
        if args.annotated_exons is None:
            msg = "Splitting the annotated transcripts into exon blocks"
            logger.info(msg)

            annotated_exons = bed_utils.split_bed12(annotated_transcripts, 
                                    num_cpus=args.num_cpus, progress_bar=True)
        else:
            msg = "Reading the annotated transcript exons"
            logger.info(msg)

            annotated_exons = bed_utils.read_bed(args.annotated_exons)

        msg = "Finding completely covered extracted ORFs"
        logger.info(msg)

        nonoverlapping_ids = bed_utils.subtract_bed(extracted_orf_exons, annotated_exons, 
                                                    min_a_overlap=1)

        m_unfiltered = extracted_orfs['id'].isin(nonoverlapping_ids)
        extracted_orfs = extracted_orfs[m_unfiltered]

        # also discard the unnecessary exons
        m_unfiltered = extracted_orf_exons['id'].isin(nonoverlapping_ids)
        extracted_orf_exons = extracted_orf_exons[m_unfiltered]

        msg = "After filtering, {} extracted ORFs remain".format(len(extracted_orfs))
        logger.info(msg)

    
    # if the nonoverlapping-label is given, annotate and remove the ORFs
    # which do not at all overlap the annotations
    if args.nonoverlapping_label is not None:
        
        nonoverlapping_ids = bed_utils.subtract_bed(
                    extracted_orfs, 
                    annotated_transcripts, 
                    exons_a=extracted_orf_exons,
                    exons_b=annotated_exons)

        m_nonoverlapping = extracted_orf_exons['id'].isin(nonoverlapping_ids)
        extracted_orf_exons = extracted_orf_exons[~m_nonoverlapping]

        m_nonoverlapping = extracted_orfs['id'].isin(nonoverlapping_ids)
        extracted_orfs.loc[m_nonoverlapping, 'orf_type'] = args.nonoverlapping_label

        msg = ("Found {} ORFs completely nonoverlapping annotated transcripts".
            format(len(nonoverlapping_ids)))
        logger.info(msg)

    msg = "Removing the annotated UTRs from the transcripts"
    logger.info(msg)
    canonical_orfs = bed_utils.retain_all_thick_only(annotated_transcripts, 
        num_cpus=args.num_cpus)

    msg = "Splitting the canonical ORFs into exons"
    logger.info(msg)
    canonical_orf_exons = bed_utils.split_bed12(canonical_orfs, 
        num_cpus=args.num_cpus, progress_bar=True)

    msg = "Extracting annotated 5' leader regions"
    logger.info(msg)
    five_prime_regions = bed_utils.retain_all_five_prime_of_thick(
        annotated_transcripts, num_cpus=args.num_cpus)

    msg = "Splitting the 5' leaders into exons"
    logger.info(msg)
    five_prime_exons = bed_utils.split_bed12(five_prime_regions, 
        num_cpus=args.num_cpus, progress_bar=True)

    msg = "Extracting annotated 3' trailer regions"
    logger.info(msg)
    three_prime_regions = bed_utils.retain_all_three_prime_of_thick(
        annotated_transcripts, num_cpus=args.num_cpus)

    msg = "Splitting the 3' trailers into exons"
    logger.info(msg)
    three_prime_exons = bed_utils.split_bed12(three_prime_regions, 
        num_cpus=args.num_cpus, progress_bar=True)

    msg = "Splitting noncoding transcripts into exons"
    logger.info(msg)

    m_no_thick_start = annotated_transcripts['thick_start'] == -1
    m_no_thick_end = annotated_transcripts['thick_end'] == -1
    m_no_thick = m_no_thick_start & m_no_thick_end
    noncoding_transcripts = annotated_transcripts[m_no_thick]

    noncoding_exons = bed_utils.split_bed12(noncoding_transcripts, 
                                    num_cpus=args.num_cpus, progress_bar=True)

    msg = "Marking canonical and extracted ORFs with the same stop codon"
    logger.info(msg)

    # first, add the true ORF end
    m_forward_canonical = canonical_orfs['strand'] == '+'
    m_reverse_canonical = canonical_orfs['strand'] == '-'

    m_forward_extracted = extracted_orfs['strand'] == '+'
    m_reverse_extracted = extracted_orfs['strand'] == '-'

    canonical_orfs['orf_end'] = canonical_orfs['end']
    canonical_orfs.loc[m_reverse_canonical, 'orf_end'] = canonical_orfs.loc[m_reverse_canonical, 'start']

    extracted_orfs['orf_end'] = extracted_orfs['end']
    extracted_orfs.loc[m_reverse_extracted, 'orf_end'] = extracted_orfs.loc[m_reverse_extracted, 'start']

    # now, find extracted ORFs with the same "orf_end" (and seqname, strand) as canonical ORFs
    merge_fields = ['seqname', 'strand', 'orf_end']
    canonical_extracted_orf_ends = canonical_orfs.merge(extracted_orfs, 
                        on=merge_fields, suffixes=['_canonical', '_extracted'])

    # now, pull this into a set
    zip_it =  zip(  canonical_extracted_orf_ends['id_canonical'], 
                    canonical_extracted_orf_ends['id_extracted'])
    canonical_extracted_matching_ends = { (c,a) for c,a in zip_it }

    msg = "Finding ORFs which exactly overlap the canonical ORFs"
    logger.info(msg)

    exact_matches = bed_utils.get_bed_overlaps(canonical_orf_exons, 
                        extracted_orf_exons, min_a_overlap=1, min_b_overlap=1)

    exact_match_orf_ids = {o.b_info for o in exact_matches}
    m_exact_orf_matches = extracted_orf_exons['id'].isin(exact_match_orf_ids)
    extracted_orf_exons = extracted_orf_exons[~m_exact_orf_matches]

    m_canonical = extracted_orfs['id'].isin(exact_match_orf_ids)
    extracted_orfs.loc[m_canonical, 'orf_type'] = 'canonical'

    msg = "Found {} canonical ORFs".format(len(exact_match_orf_ids))
    logger.info(msg)

    msg = "Finding ORFs which are extended versions of the canonical ORFs"
    logger.info(msg)

    extended_matches = bed_utils.get_bed_overlaps(canonical_orf_exons, 
                                        extracted_orf_exons, min_a_overlap=1)

    # make sure the "end"s match before calling something an extended match
    extended_match_ids = {m.b_info for m in tqdm.tqdm(extended_matches) 
                            if (m.a_info, m.b_info) in canonical_extracted_matching_ends}

    m_extended_matches = extracted_orf_exons['id'].isin(extended_match_ids)
    extracted_orf_exons = extracted_orf_exons[~m_extended_matches]

    m_canonical_extended = extracted_orfs['id'].isin(extended_match_ids)

    l = "{}canonical_extended".format(args.label_prefix)
    extracted_orfs.loc[m_canonical_extended, 'orf_type'] = l

    msg = "Found {} canonical_extended ORFs".format(len(extended_match_ids))
    logger.info(msg)

    msg = "Finding ORFs which are truncated versions of the canonical ORFs"
    logger.info(msg)

    truncated_matches = bed_utils.get_bed_overlaps(canonical_orf_exons, 
                                        extracted_orf_exons, min_b_overlap=1)

    # make sure the "end"s match before calling something a truncated match
    truncated_match_ids = {m.b_info for m in tqdm.tqdm(truncated_matches) 
                            if (m.a_info, m.b_info) in canonical_extracted_matching_ends}

    m_truncated_matches = extracted_orf_exons['id'].isin(truncated_match_ids)
    extracted_orf_exons = extracted_orf_exons[~m_truncated_matches]

    m_canonical_truncated = extracted_orfs['id'].isin(truncated_match_ids)
    
    l = "{}canonical_truncated".format(args.label_prefix)
    extracted_orfs.loc[m_canonical_truncated, 'orf_type'] = l

    msg = "Found {} canonical_truncated ORFs".format(len(truncated_match_ids))
    logger.info(msg)

    msg = ("Labeling ORFs which are completely covered by a canonical ORF but "
                    "do not share its stop codon")
    logger.info(msg)

    # anything in "truncated matches" which *does not* share a stop codon with 
    # the match is a "within" orf
    within_ids = {m.b_info for m in truncated_matches if m.b_info not in truncated_match_ids}

    m_within_matches = extracted_orf_exons['id'].isin(within_ids)
    extracted_orf_exons = extracted_orf_exons[~m_within_matches]

    m_within = extracted_orfs['id'].isin(within_ids)
    
    l = "{}within".format(args.label_prefix)
    extracted_orfs.loc[m_within, 'orf_type'] = l

    msg = "Found {} within ORFs".format(len(within_ids))
    logger.info(msg)
        
    msg = "Finding out-of-frame overlaps"
    logger.info(msg)
    out_of_frame_matches = bed_utils.get_bed_overlaps(canonical_orf_exons, 
                                                        extracted_orf_exons)

    msg = "Finding leader overlaps"
    logger.info(msg)

    leader_matches = bed_utils.get_bed_overlaps(five_prime_exons, 
                                                    extracted_orf_exons)

    msg = "Finding trailer overlaps"
    logger.info(msg)

    trailer_matches = bed_utils.get_bed_overlaps(three_prime_exons, 
                                                    extracted_orf_exons)

    msg = ("Labeling ORFs which have (out-of-frame) overlaps with both a "
        "canonical ORF and annotated leaders or trailers")
    logger.info(msg)

    out_of_frame_ids = {m.b_info for m in out_of_frame_matches}
    leader_ids = {m.b_info for m in leader_matches}
    trailer_ids = {m.b_info for m in trailer_matches}

    leader_overlap_ids = out_of_frame_ids & leader_ids
    trailer_overlap_ids = out_of_frame_ids & trailer_ids

    m_leader_overlap_matches = extracted_orf_exons['id'].isin(leader_overlap_ids)
    extracted_orf_exons = extracted_orf_exons[~m_leader_overlap_matches]

    m_trailer_overlap_matches = extracted_orf_exons['id'].isin(trailer_overlap_ids)
    extracted_orf_exons = extracted_orf_exons[~m_trailer_overlap_matches]

    m_five_prime_overlap = extracted_orfs['id'].isin(leader_overlap_ids)
    
    l = "{}five_prime_overlap".format(args.label_prefix)
    extracted_orfs.loc[m_five_prime_overlap, 'orf_type'] = l

    m_three_prime_overlap = extracted_orfs['id'].isin(trailer_overlap_ids)
    
    l = "{}three_prime_overlap".format(args.label_prefix)
    extracted_orfs.loc[m_three_prime_overlap, 'orf_type'] = l

    msg = "Found {} five_prime_overlap ORFs".format(len(leader_overlap_ids))
    logger.info(msg)

    msg = "Found {} three_prime_overlap ORFs".format(len(trailer_overlap_ids))
    logger.info(msg)

    msg = "Finding ORFs completely within 5' leaders"
    logger.info(msg)

    leader_matches = bed_utils.get_bed_overlaps(five_prime_exons, 
                                        extracted_orf_exons, min_b_overlap=1)
    leader_ids = {m.b_info for m in leader_matches}

    m_leader_matches = extracted_orf_exons['id'].isin(leader_ids)
    extracted_orf_exons = extracted_orf_exons[~m_leader_matches]

    m_five_prime = extracted_orfs['id'].isin(leader_ids)
    
    l = "{}five_prime".format(args.label_prefix)
    extracted_orfs.loc[m_five_prime, 'orf_type'] = l

    msg = "Found {} five_prime ORFs".format(len(leader_ids))
    logger.info(msg)

    msg = "Finding ORFs completely within 3' trailers"
    logger.info(msg)

    trailer_matches = bed_utils.get_bed_overlaps(three_prime_exons, 
                                        extracted_orf_exons, min_b_overlap=1)
    trailer_ids = {m.b_info for m in trailer_matches}

    m_trailer_matches = extracted_orf_exons['id'].isin(trailer_ids)
    extracted_orf_exons = extracted_orf_exons[~m_trailer_matches]

    m_three_prime = extracted_orfs['id'].isin(trailer_ids)
    
    l = "{}three_prime".format(args.label_prefix)
    extracted_orfs.loc[m_three_prime, 'orf_type'] = l

    msg = "Found {} three_prime ORFs".format(len(trailer_ids))
    logger.info(msg)

    msg = "Finding ORFs completely within annotated, noncoding transcripts"
    logger.info(msg)

    noncoding_matches = bed_utils.get_bed_overlaps(noncoding_exons, 
                                        extracted_orf_exons, min_b_overlap=1)

    noncoding_ids = {m.b_info for m in noncoding_matches}

    m_noncoding_matches = extracted_orf_exons['id'].isin(noncoding_ids)
    extracted_orf_exons = extracted_orf_exons[~m_noncoding_matches]

    m_noncoding = extracted_orfs['id'].isin(noncoding_ids)
    
    l = "{}noncoding".format(args.label_prefix)
    extracted_orfs.loc[m_noncoding, 'orf_type'] = l

    msg = "Found {} noncoding ORFs".format(len(noncoding_ids))
    logger.info(msg)

    # all of the remaining ORFs fall into the "suspect" category
    suspect_ids = {orf_id for orf_id in extracted_orf_exons['id']}

    m_suspect = extracted_orfs['id'].isin(suspect_ids)
    
    l = "{}suspect".format(args.label_prefix)
    extracted_orfs.loc[m_suspect, 'orf_type'] = l

    msg = "Found {} \"suspect\" ORFs".format(len(suspect_ids))
    logger.info(msg)

    m_no_orf_type = extracted_orfs['orf_type'].isnull()

    msg = "Found {} unlabeled ORFs".format(sum(m_no_orf_type))
    logger.info(msg)

    msg = "Writing ORFs with types to disk"
    logger.info(msg)

    fields = bed_utils.bed12_field_names + ['orf_num', 'orf_len', 'orf_type']
    extracted_orfs = extracted_orfs[fields]
    extracted_orfs = bed_utils.sort(extracted_orfs)

    bed_utils.write_bed(extracted_orfs, args.out)

if __name__ == '__main__':
    main()
