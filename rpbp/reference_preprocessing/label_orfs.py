#! /usr/bin/env python3

import argparse
import logging
import tqdm

import pandas as pd

import misc.logging_utils as logging_utils
import misc.parallel as parallel

import bio_utils.bed_utils as bed_utils

logger = logging.getLogger(__name__)

default_num_cpus = 1

default_annotated_exons = None
default_label_prefix = ""
default_nonoverlapping_label = None


def change_orf_id(row, match_pairs):
    """Make sure the label matches the transcript part of the
    orf id, i.e. we want the label to be consistent with the
    transcript to which the ORF is associated."""

    orf_id = row.id
    if row.is_duplicated:
        split_id = orf_id.rpartition('_')
        transcript_id = split_id[0]
        if (transcript_id, orf_id) not in match_pairs:
            # just pick the first one, all associated transcripts
            # share the same relationship with this ORF
            transcript_id = [pair[0] for pair in match_pairs if pair[1] == orf_id][0]
            orf_id = transcript_id + ''.join(split_id[1:])

    return orf_id


def check_orf_ids(orf_group, type_to_matches):

    orf_type = orf_group['orf_type'].iloc[0]
    match_pairs = type_to_matches[orf_type]

    orf_group['id'] = orf_group.apply(change_orf_id,
                                      match_pairs=match_pairs,
                                      axis=1)

    return orf_group


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='''Label the ORFs based on their exon structure 
                                     wrt the annotated transcripts. The extracted_orfs input 
                                     file can be overwritten.''')
    parser.add_argument('annotated_transcripts', help='''The annotated transcripts for the genome
        in BED12+ format.''')
    parser.add_argument('extracted_orfs', help='''The ORFs extracted from the transcripts 
        in BED12+ format.''')
    parser.add_argument('orf_exons', help='''The exon blocks for the ORFs, in BED6+ format, obtained
        from "split-bed12-blocks".''')
    parser.add_argument('out', help='''The output (BED12+.gz) file.''')
    parser.add_argument('-p', '--num-cpus', help='''The number of CPUs to use to perform
        BED operations.''', type=int, default=default_num_cpus)
    parser.add_argument('-f', '--filter', help='''If this flag is given, then ORFs
        which are completely covered by an annotated transcript are discarded. Use to filter 
        uninteresting ORFs from a de novo assembly.''', action='store_true')
    parser.add_argument('-e', '--annotated-exons', help='''The annotated transcript 
        exons can be passed with this option. If they are not given, they will be 
        split from the annotated transcripts.''', default=default_annotated_exons)
    parser.add_argument('-n', '--nonoverlapping-label', help='''If this option is given, 
        then the ORFs which do not overlap the annotated transcripts at all will be given this label.
        By default, remaining oof overlapping ORFs are assigned the "overlap" label.
        If not given, the ORFs outside of annotated regions are labeled as "suspect".''',
                        default=default_nonoverlapping_label)
    parser.add_argument('-l', '--label-prefix', help='''This string is prepended to all labels 
        assigned to ORFs, e.g. to indicate ORFs from a de novo assembly (Rp-Bp assigns the label
        "novel" to these, however the string is not prepended to "canonical ORFs").''',
                        default=default_label_prefix)
    parser.add_argument('-s', '--skip-check', help='''If this flag is given, then the 
        transcript part of the orf_id is not checked for consistency with the label. Use
        when the orf_id is not associated with a transcript_id, such as in a de novo assembly.
        If given, then the ORF exons file will be re-written accordingly.''', action='store_true')

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

    # get the annotated transcript exons
    if args.annotated_exons is None:
        msg = "Splitting the annotated transcripts into exon blocks"
        logger.info(msg)

        annotated_exons = bed_utils.split_bed12(annotated_transcripts,
                                                num_cpus=args.num_cpus,
                                                progress_bar=True)
    else:
        msg = "Reading the annotated transcript exons"
        logger.info(msg)

        annotated_exons = bed_utils.read_bed(args.annotated_exons)

    # filter out the ORFs that are entirely within annotated transcripts
    if args.filter:
        msg = ("Removing extracted ORFs which are completely covered by the "
               "annotated transcripts")
        logger.info(msg)
        msg = "Finding completely covered extracted ORFs"
        logger.info(msg)

        nonoverlapping_ids = bed_utils.subtract_bed(extracted_orf_exons,
                                                    annotated_exons,
                                                    min_a_overlap=1)
        m_unfiltered = extracted_orfs['id'].isin(nonoverlapping_ids)
        extracted_orfs = extracted_orfs[m_unfiltered]
        # discard the unnecessary exons
        m_unfiltered = extracted_orf_exons['id'].isin(nonoverlapping_ids)
        extracted_orf_exons = extracted_orf_exons[m_unfiltered]

        msg = "After filtering, {} extracted ORFs remain".format(len(extracted_orfs))
        logger.info(msg)

    # annotate and remove the ORFs which do not at all overlap with the annotations
    if args.nonoverlapping_label is not None:
        nonoverlapping_ids = bed_utils.subtract_bed(extracted_orfs,
                                                    annotated_transcripts,
                                                    exons_a=extracted_orf_exons,
                                                    exons_b=annotated_exons)
        m_nonoverlapping = extracted_orf_exons['id'].isin(nonoverlapping_ids)
        extracted_orf_exons = extracted_orf_exons[~m_nonoverlapping]
        m_nonoverlapping = extracted_orfs['id'].isin(nonoverlapping_ids)
        extracted_orfs.loc[m_nonoverlapping, 'orf_type'] = args.nonoverlapping_label

        msg = ("Found {} ORFs completely non-overlapping annotated transcripts".
               format(len(nonoverlapping_ids)))
        logger.info(msg)

    msg = "Removing the annotated UTRs from the transcripts"
    logger.info(msg)
    canonical_orfs = bed_utils.retain_all_thick_only(annotated_transcripts,
                                                     num_cpus=args.num_cpus)

    msg = "Splitting the canonical ORFs into exons"
    logger.info(msg)
    canonical_orf_exons = bed_utils.split_bed12(canonical_orfs,
                                                num_cpus=args.num_cpus,
                                                progress_bar=True)

    msg = "Extracting annotated 5' leader regions"
    logger.info(msg)
    five_prime_regions = bed_utils.retain_all_five_prime_of_thick(
        annotated_transcripts, num_cpus=args.num_cpus)

    if len(five_prime_regions) == 0:
        msg = "No annotated 5' leader regions were found"
        logger.warning(msg)

    msg = "Splitting the 5' leaders into exons"
    logger.info(msg)
    five_prime_exons = bed_utils.split_bed12(five_prime_regions,
                                             num_cpus=args.num_cpus,
                                             progress_bar=True)

    msg = "Extracting annotated 3' trailer regions"
    logger.info(msg)
    three_prime_regions = bed_utils.retain_all_three_prime_of_thick(
        annotated_transcripts, num_cpus=args.num_cpus)

    if len(three_prime_regions) == 0:
        msg = "No annotated 3' trailer regions were found"
        logger.warning(msg)

    msg = "Splitting the 3' trailers into exons"
    logger.info(msg)
    three_prime_exons = bed_utils.split_bed12(three_prime_regions,
                                              num_cpus=args.num_cpus,
                                              progress_bar=True)

    msg = "Splitting non-coding transcripts into exons"
    logger.info(msg)

    m_no_thick_start = annotated_transcripts['thick_start'] == -1
    m_no_thick_end = annotated_transcripts['thick_end'] == -1
    m_no_thick = m_no_thick_start & m_no_thick_end
    noncoding_transcripts = annotated_transcripts[m_no_thick]
    noncoding_exons = bed_utils.split_bed12(noncoding_transcripts,
                                            num_cpus=args.num_cpus,
                                            progress_bar=True)

    # First, remove all in-frame (canonical, canonical variants), and also within (oof) ORFs

    msg = "Marking canonical and extracted ORFs with the same stop codon"
    logger.info(msg)

    # first, add the "true" ORF end
    m_reverse_canonical = canonical_orfs['strand'] == '-'
    canonical_orfs['orf_end'] = canonical_orfs['end']
    canonical_orfs.loc[m_reverse_canonical, 'orf_end'] = canonical_orfs.loc[m_reverse_canonical, 'start']

    m_reverse_extracted = extracted_orfs['strand'] == '-'
    extracted_orfs['orf_end'] = extracted_orfs['end']
    extracted_orfs.loc[m_reverse_extracted, 'orf_end'] = extracted_orfs.loc[m_reverse_extracted, 'start']

    # then, find extracted ORFs with the same "orf_end" (and seqname, strand) as canonical ORFs
    merge_fields = ['seqname', 'strand', 'orf_end']
    canonical_extracted_orf_ends = canonical_orfs.merge(extracted_orfs,
                                                        on=merge_fields,
                                                        suffixes=['_canonical', '_extracted'])

    # finally, pull this into a set
    zip_it = zip(canonical_extracted_orf_ends['id_canonical'],
                 canonical_extracted_orf_ends['id_extracted'])
    canonical_extracted_matching_ends = {(c, a) for c, a in zip_it}

    msg = "Finding ORFs which exactly overlap the canonical ORFs"
    logger.info(msg)

    exact_matches = bed_utils.get_bed_overlaps(canonical_orf_exons,
                                               extracted_orf_exons,
                                               min_a_overlap=1,
                                               min_b_overlap=1)

    exact_match_orf_ids = {m.b_info for m in exact_matches}
    m_exact_orf_matches = extracted_orf_exons['id'].isin(exact_match_orf_ids)
    extracted_orf_exons = extracted_orf_exons[~m_exact_orf_matches]
    m_canonical = extracted_orfs['id'].isin(exact_match_orf_ids)
    extracted_orfs.loc[m_canonical, 'orf_type'] = 'canonical'

    msg = "Found {} canonical ORFs".format(len(exact_match_orf_ids))
    logger.info(msg)

    msg = "Finding ORFs which are extended versions of the canonical ORFs"
    logger.info(msg)

    extended_matches = bed_utils.get_bed_overlaps(canonical_orf_exons,
                                                  extracted_orf_exons,
                                                  min_a_overlap=1)

    # make sure the "end"s match before calling something an extended match
    extended_match_pairs = {(m.a_info, m.b_info) for m in tqdm.tqdm(extended_matches)
                            if (m.a_info, m.b_info) in canonical_extracted_matching_ends}
    extended_match_ids = {pair[1] for pair in extended_match_pairs}

    m_extended_matches = extracted_orf_exons['id'].isin(extended_match_ids)
    extracted_orf_exons = extracted_orf_exons[~m_extended_matches]
    m_canonical_extended = extracted_orfs['id'].isin(extended_match_ids)

    label = "{}canonical_extended".format(args.label_prefix)
    extracted_orfs.loc[m_canonical_extended, 'orf_type'] = label

    type_to_matches = {label: extended_match_pairs}

    msg = "Found {} canonical_extended ORFs".format(len(extended_match_ids))
    logger.info(msg)

    msg = "Finding ORFs which are truncated versions of the canonical ORFs"
    logger.info(msg)

    truncated_matches = bed_utils.get_bed_overlaps(canonical_orf_exons,
                                                   extracted_orf_exons,
                                                   min_b_overlap=1)

    # make sure the "end"s match before calling something a truncated match
    truncated_match_pairs = {(m.a_info, m.b_info) for m in tqdm.tqdm(truncated_matches)
                             if (m.a_info, m.b_info) in canonical_extracted_matching_ends}
    truncated_match_ids = {pair[1] for pair in truncated_match_pairs}

    m_truncated_matches = extracted_orf_exons['id'].isin(truncated_match_ids)
    extracted_orf_exons = extracted_orf_exons[~m_truncated_matches]
    m_canonical_truncated = extracted_orfs['id'].isin(truncated_match_ids)

    label = "{}canonical_truncated".format(args.label_prefix)
    extracted_orfs.loc[m_canonical_truncated, 'orf_type'] = label

    type_to_matches[label] = truncated_match_pairs

    msg = "Found {} canonical_truncated ORFs".format(len(truncated_match_ids))
    logger.info(msg)

    msg = ("Labeling ORFs which are completely covered by a canonical ORF but "
           "are out-of-frame")
    logger.info(msg)

    # anything in "truncated matches" which *does not* share a stop codon with
    # the match is a "within" ORF ... not necessarily, see
    # https://github.com/dieterich-lab/rp-bp/issues/96
    within_pairs = {(m.a_info, m.b_info) for m in truncated_matches
                    if m.b_info not in truncated_match_ids}
    within_ids = {pair[1] for pair in within_pairs}

    m_within_matches = extracted_orf_exons['id'].isin(within_ids)
    extracted_orf_exons = extracted_orf_exons[~m_within_matches]
    m_within = extracted_orfs['id'].isin(within_ids)

    if args.nonoverlapping_label is None:
        label = "{}within".format(args.label_prefix)
        type_to_matches[label] = within_pairs

        msg = "Found {} within ORFs".format(len(within_ids))
        logger.info(msg)
    else:
        label = "{}overlap".format(args.label_prefix)

    extracted_orfs.loc[m_within, 'orf_type'] = label

    # then, find ORFs that are overlapping, or that are entirely within the UTRs

    msg = "Finding all overlap matches"
    logger.info(msg)
    out_of_frame_matches = bed_utils.get_bed_overlaps(canonical_orf_exons,
                                                      extracted_orf_exons)

    leader_matches = bed_utils.get_bed_overlaps(five_prime_exons,
                                                extracted_orf_exons)

    trailer_matches = bed_utils.get_bed_overlaps(three_prime_exons,
                                                 extracted_orf_exons)

    msg = ("Labeling ORFs which have (out-of-frame) overlaps with both a "
           "canonical ORF and annotated leaders or trailers")
    logger.info(msg)

    # make sure up-/downstream overlaps are unique,
    # in cases where an ORF overlaps both the 5'UTR and
    # the 3'UTR of different same sense overlapping transcripts,
    # it is assigned by default to the downstream overlap
    # for de novo, everything is labeled as overlap

    leader_match_pairs = {(m.a_info, m.b_info) for m in leader_matches}
    trailer_match_pairs = {(m.a_info, m.b_info) for m in trailer_matches}

    if args.nonoverlapping_label is None:

        # for standard assembly, we also need to make sure
        # all overlap matches are fully contained within the
        # transcript structure
        transcript_matches = bed_utils.get_bed_overlaps(annotated_exons,
                                                        extracted_orf_exons,
                                                        min_b_overlap=1)

        transcript_match_pairs = {(m.a_info, m.b_info) for m in transcript_matches}

        leader_overlap_pairs = {(m.a_info, m.b_info) for m in out_of_frame_matches
                                if (m.a_info, m.b_info) in leader_match_pairs
                                and (m.a_info, m.b_info) not in trailer_match_pairs
                                and (m.a_info, m.b_info) in transcript_match_pairs}

        trailer_overlap_pairs = {(m.a_info, m.b_info) for m in out_of_frame_matches
                                 if (m.a_info, m.b_info) in trailer_match_pairs
                                 and (m.a_info, m.b_info) not in leader_match_pairs
                                 and (m.a_info, m.b_info) in transcript_match_pairs}

        # we do not assign preference in cases where the ORF overlaps both sides
        # of the coding sequence on the same transcript, any ORF
        # satisfying both will be labeled simply as overlap
        overlap_pairs = {(m.a_info, m.b_info) for m in out_of_frame_matches
                         if (m.a_info, m.b_info) in leader_match_pairs
                         and (m.a_info, m.b_info) in trailer_match_pairs
                         and (m.a_info, m.b_info) in transcript_match_pairs}

        overlap_ids = {pair[1] for pair in overlap_pairs}

        trailer_overlap_pairs = {pair for pair in trailer_overlap_pairs
                                 if pair[1] not in overlap_ids}
        trailer_overlap_ids = {pair[1] for pair in trailer_overlap_pairs}

        leader_overlap_pairs = {pair for pair in leader_overlap_pairs
                                if pair[1] not in trailer_overlap_ids
                                and pair[1] not in overlap_ids}
        leader_overlap_ids = {pair[1] for pair in leader_overlap_pairs}

        m_overlap_matches = extracted_orf_exons['id'].isin(overlap_ids)
        extracted_orf_exons = extracted_orf_exons[~m_overlap_matches]

        m_leader_overlap_matches = extracted_orf_exons['id'].isin(leader_overlap_ids)
        extracted_orf_exons = extracted_orf_exons[~m_leader_overlap_matches]

        m_trailer_overlap_matches = extracted_orf_exons['id'].isin(trailer_overlap_ids)
        extracted_orf_exons = extracted_orf_exons[~m_trailer_overlap_matches]

        m_five_prime_overlap = extracted_orfs['id'].isin(leader_overlap_ids)
        label = "{}five_prime_overlap".format(args.label_prefix)
        extracted_orfs.loc[m_five_prime_overlap, 'orf_type'] = label

        type_to_matches[label] = leader_overlap_pairs

        m_three_prime_overlap = extracted_orfs['id'].isin(trailer_overlap_ids)
        label = "{}three_prime_overlap".format(args.label_prefix)
        extracted_orfs.loc[m_three_prime_overlap, 'orf_type'] = label

        type_to_matches[label] = trailer_overlap_pairs

        msg = "Found {} five_prime_overlap ORFs".format(len(leader_overlap_ids))
        logger.info(msg)
        msg = "Found {} three_prime_overlap ORFs".format(len(trailer_overlap_ids))
        logger.info(msg)

    else:

        overlap_pairs = {(m.a_info, m.b_info) for m in out_of_frame_matches}
        overlap_pairs |= leader_match_pairs
        overlap_pairs |= trailer_match_pairs
        overlap_pairs |= within_pairs
        overlap_ids = {pair[1] for pair in overlap_pairs}

        m_overlap_matches = extracted_orf_exons['id'].isin(overlap_ids)
        extracted_orf_exons = extracted_orf_exons[~m_overlap_matches]

    m_overlap = extracted_orfs['id'].isin(overlap_ids)
    label = "{}overlap".format(args.label_prefix)
    extracted_orfs.loc[m_overlap, 'orf_type'] = label

    type_to_matches[label] = overlap_pairs

    msg = "Found {} overlap ORFs".format(len(overlap_ids))
    logger.info(msg)

    msg = "Finding ORFs completely within 5' or 3' leaders"
    logger.info(msg)

    # here, upstream ORFs are labeled first and immediately removed,
    # so there is no double match, we just need to make sure below
    # that the transcript part of the "orf_id" is consistent with the label.

    leader_matches = bed_utils.get_bed_overlaps(five_prime_exons,
                                                extracted_orf_exons,
                                                min_b_overlap=1)

    leader_pairs = {(m.a_info, m.b_info) for m in leader_matches}
    leader_ids = {pair[1] for pair in leader_pairs}

    m_leader_matches = extracted_orf_exons['id'].isin(leader_ids)
    extracted_orf_exons = extracted_orf_exons[~m_leader_matches]

    m_five_prime = extracted_orfs['id'].isin(leader_ids)
    label = "{}five_prime".format(args.label_prefix)
    extracted_orfs.loc[m_five_prime, 'orf_type'] = label

    type_to_matches[label] = leader_pairs

    msg = "Found {} five_prime ORFs".format(len(leader_ids))
    logger.info(msg)

    trailer_matches = bed_utils.get_bed_overlaps(three_prime_exons,
                                                 extracted_orf_exons,
                                                 min_b_overlap=1)

    trailer_pairs = {(m.a_info, m.b_info) for m in trailer_matches}
    trailer_ids = {pair[1] for pair in trailer_pairs}

    m_trailer_matches = extracted_orf_exons['id'].isin(trailer_ids)
    extracted_orf_exons = extracted_orf_exons[~m_trailer_matches]

    m_three_prime = extracted_orfs['id'].isin(trailer_ids)
    label = "{}three_prime".format(args.label_prefix)
    extracted_orfs.loc[m_three_prime, 'orf_type'] = label

    type_to_matches[label] = trailer_pairs

    msg = "Found {} three_prime ORFs".format(len(trailer_ids))
    logger.info(msg)

    msg = "Finding ORFs completely within annotated, noncoding transcripts"
    logger.info(msg)

    noncoding_matches = bed_utils.get_bed_overlaps(noncoding_exons,
                                                   extracted_orf_exons,
                                                   min_b_overlap=1)

    noncoding_pairs = {(m.a_info, m.b_info) for m in noncoding_matches}
    noncoding_ids = {pair[1] for pair in noncoding_pairs}

    m_noncoding_matches = extracted_orf_exons['id'].isin(noncoding_ids)
    extracted_orf_exons = extracted_orf_exons[~m_noncoding_matches]

    m_noncoding = extracted_orfs['id'].isin(noncoding_ids)
    label = "{}noncoding".format(args.label_prefix)
    extracted_orfs.loc[m_noncoding, 'orf_type'] = label

    type_to_matches[label] = noncoding_pairs

    msg = "Found {} noncoding ORFs".format(len(noncoding_ids))
    logger.info(msg)

    # all of the remaining ORFs fall into the "suspect" category
    suspect_ids = {orf_id for orf_id in extracted_orf_exons['id']}

    m_suspect = extracted_orfs['id'].isin(suspect_ids)
    label = "{}suspect".format(args.label_prefix)
    extracted_orfs.loc[m_suspect, 'orf_type'] = label

    msg = "Remaining {} ORFs labeled as suspect".format(len(suspect_ids))
    logger.info(msg)

    m_no_orf_type = extracted_orfs['orf_type'].isnull()
    msg = "Found {} unlabeled ORFs".format(sum(m_no_orf_type))
    logger.info(msg)

    msg = "Now checking all labels against orf_ids"
    logger.info(msg)

    if not args.skip_check:

        grouped = extracted_orfs.groupby('orf_type')
        grouped_canonical = grouped.get_group('canonical')
        orfs_group = extracted_orfs.drop(grouped_canonical.index).groupby('orf_type')

        relabeled_orfs = parallel.apply_parallel_groups(orfs_group, args.num_cpus,
                                                        check_orf_ids, type_to_matches,
                                                        progress_bar=True)
        relabeled_orfs.append(grouped_canonical)
        extracted_orfs = pd.concat(relabeled_orfs)

    msg = "Writing ORFs with types to disk"
    logger.info(msg)

    extracted_orfs = bed_utils.sort(extracted_orfs)
    additional_columns = ['orf_num', 'orf_len', 'orf_type']
    if 'assoc_trx' in extracted_orfs.columns:
        additional_columns.extend(['assoc_trx'])
    fields = bed_utils.bed12_field_names + additional_columns
    extracted_orfs = extracted_orfs[fields]

    bed_utils.write_bed(extracted_orfs, args.out)


if __name__ == '__main__':
    main()
