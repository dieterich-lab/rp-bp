#! /usr/bin/env python3

"""This script labels the ORFs based on their exon
transcript structure with respect to annotated coding sequences
"""

import argparse
import logging

import pbiotools.misc.logging_utils as logging_utils
import pbiotools.utils.bed_utils as bed_utils

from rpbp.defaults import default_num_cpus, orf_type_name_map

logger = logging.getLogger(__name__)


def get_labels(args):
    # standard labels
    labels = orf_type_name_map.copy()
    # standalone: outside the Rp-Bp pipeline
    # we need to make sure [--label-prefix] and [--nonoverlapping-label]
    # are contained in the keys
    # NOTE: This is merely to avoid a KeyError, as these options are designed to handle
    # ORFs from a de novo assembly, and not simply as a general label tag.
    if (
        args.nonoverlapping_label is not None
        and not args.nonoverlapping_label == "novel"
    ):
        labels[args.nonoverlapping_label] = args.nonoverlapping_label
    if len(args.label_prefix) > 0 and not args.label_prefix == "novel_":
        for key, val in orf_type_name_map.items():
            # replace standard labels
            lkey = key.replace("novel_", args.label_prefix)
            labels[lkey] = lkey
    return labels


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Label the ORFs based on their transcript
        exon structure wrt the annotated transcripts.""",
    )

    parser.add_argument(
        "annotated_transcripts",
        help="""The annotated transcripts for the genome
        in BED12+ format.""",
    )

    parser.add_argument(
        "extracted_orfs",
        help="""The ORFs extracted from the transcripts
        in BED12+ format.""",
    )

    parser.add_argument("out", help="""The output (BED12+.gz) file.""")

    parser.add_argument(
        "-e",
        "--annotated-exons",
        help="""The annotated transcript
        exons can be passed with this option. If they are not given, they will be
        split from the annotated transcripts.""",
        default=None,
    )

    parser.add_argument(
        "-o",
        "--orf-exons",
        help="""The exon blocks for the ORFs, in BED6+ format,
        obtained from "split-bed12-blocks". If they are not given, they will be split from the
        extracted ORFs.""",
        default=None,
    )

    parser.add_argument(
        "-n",
        "--nonoverlapping-label",
        help="""If this option is given,
        then the ORFs which do not overlap the annotated transcripts at all will be given this label.
        By default, remaining oof overlapping ORFs are assigned the "overlap" label.
        If not given, the ORFs outside of annotated regions are labeled as "suspect".""",
        default=None,
    )

    parser.add_argument(
        "-l",
        "--label-prefix",
        help="""This string is prepended to all labels
        assigned to ORFs, e.g. to indicate ORFs from a de novo assembly (Rp-Bp assigns the label
        "novel" to these, however the string is not prepended to "canonical ORFs").""",
        default="",
    )

    parser.add_argument(
        "-f",
        "--filter",
        help="""If this flag is given, then ORFs
        which are completely covered by an annotated transcript are discarded. Use to filter
        uninteresting ORFs from a de novo assembly.""",
        action="store_true",
    )

    parser.add_argument(
        "-p",
        "--num-cpus",
        help="""The number of CPUs to use to perform
            BED operations.""",
        type=int,
        default=default_num_cpus,
    )

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "Getting the labels"
    logger.info(msg)

    labels = get_labels(args)

    msg = "Reading annotated transcripts"
    logger.info(msg)
    annotated_transcripts = bed_utils.read_bed(args.annotated_transcripts)

    # get the annotated transcript exons
    if args.annotated_exons is None:
        msg = "Splitting the annotated transcripts into exon blocks"
        logger.info(msg)

        annotated_exons = bed_utils.split_bed12(
            annotated_transcripts, num_cpus=args.num_cpus, progress_bar=True
        )
    else:
        msg = "Reading the annotated transcript exons"
        logger.info(msg)

        annotated_exons = bed_utils.read_bed(args.annotated_exons)

    msg = "Reading extracted ORFs"
    logger.info(msg)
    extracted_orfs = bed_utils.read_bed(args.extracted_orfs)

    if args.orf_exons is None:
        msg = "Splitting the extracted ORFs into exon blocks"
        logger.info(msg)
        extracted_orf_exons = bed_utils.split_bed12(
            extracted_orfs, num_cpus=args.num_cpus, progress_bar=True
        )
    else:
        msg = "Reading the extracted ORFs exons"
        logger.info(msg)
        extracted_orf_exons = bed_utils.read_bed(args.orf_exons)

    msg = "Found {} extracted ORFs with {} exons".format(
        len(extracted_orfs), len(extracted_orf_exons)
    )
    logger.debug(msg)

    # filter out the ORFs that are entirely within annotated transcripts
    if args.filter:
        msg = "Removing ORFs which are completely covered by the annotated transcripts"
        logger.info(msg)

        nonoverlapping_ids = bed_utils.subtract_bed(
            extracted_orf_exons, annotated_exons, min_a_overlap=1
        )
        m_unfiltered = extracted_orfs["id"].isin(nonoverlapping_ids)
        extracted_orfs = extracted_orfs[m_unfiltered]
        # discard the unnecessary exons
        m_unfiltered = extracted_orf_exons["id"].isin(nonoverlapping_ids)
        extracted_orf_exons = extracted_orf_exons[m_unfiltered]

        msg = "After filtering, {} extracted ORFs remain".format(len(extracted_orfs))
        logger.info(msg)

    # annotate and remove the ORFs which do not at all overlap the annotations
    if args.nonoverlapping_label is not None:
        nonoverlapping_ids = bed_utils.subtract_bed(
            extracted_orfs,
            annotated_transcripts,
            exons_a=extracted_orf_exons,
            exons_b=annotated_exons,
        )
        m_nonoverlapping = extracted_orf_exons["id"].isin(nonoverlapping_ids)
        extracted_orf_exons = extracted_orf_exons[~m_nonoverlapping]
        m_nonoverlapping = extracted_orfs["id"].isin(nonoverlapping_ids)
        extracted_orfs.loc[m_nonoverlapping, "orf_type"] = labels[
            args.nonoverlapping_label
        ]

        msg = "Found {} ORFs completely non-overlapping annotated transcripts".format(
            len(nonoverlapping_ids)
        )
        logger.info(msg)

    msg = "Removing the annotated UTRs from the transcripts"
    logger.info(msg)
    canonical_orfs = bed_utils.retain_all_thick_only(
        annotated_transcripts, num_cpus=args.num_cpus
    )

    msg = "Splitting the canonical ORFs into exons"
    logger.info(msg)
    canonical_orf_exons = bed_utils.split_bed12(
        canonical_orfs, num_cpus=args.num_cpus, progress_bar=True
    )

    msg = "Extracting annotated 5' leader regions"
    logger.info(msg)
    five_prime_regions = bed_utils.retain_all_five_prime_of_thick(
        annotated_transcripts, num_cpus=args.num_cpus
    )

    if len(five_prime_regions) == 0:
        msg = "No annotated 5' leader regions were found"
        logger.warning(msg)

    msg = "Splitting the 5' leaders into exons"
    logger.info(msg)
    five_prime_exons = bed_utils.split_bed12(
        five_prime_regions, num_cpus=args.num_cpus, progress_bar=True
    )

    msg = "Extracting annotated 3' trailer regions"
    logger.info(msg)
    three_prime_regions = bed_utils.retain_all_three_prime_of_thick(
        annotated_transcripts, num_cpus=args.num_cpus
    )

    if len(three_prime_regions) == 0:
        msg = "No annotated 3' trailer regions were found"
        logger.warning(msg)

    msg = "Splitting the 3' trailers into exons"
    logger.info(msg)
    three_prime_exons = bed_utils.split_bed12(
        three_prime_regions, num_cpus=args.num_cpus, progress_bar=True
    )

    msg = "Splitting non-coding transcripts into exons"
    logger.info(msg)

    m_no_thick_start = annotated_transcripts["thick_start"] == -1
    m_no_thick_end = annotated_transcripts["thick_end"] == -1
    m_no_thick = m_no_thick_start & m_no_thick_end
    noncoding_transcripts = annotated_transcripts[m_no_thick]
    noncoding_exons = bed_utils.split_bed12(
        noncoding_transcripts, num_cpus=args.num_cpus, progress_bar=True
    )

    # First, remove all in-frame (canonical, canonical variants), and also within and oof ORFs

    msg = "Marking canonical and extracted ORFs with the same stop codon"
    logger.info(msg)

    # first, add the "true" ORF end
    m_reverse_canonical = canonical_orfs["strand"] == "-"
    canonical_orfs["orf_end"] = canonical_orfs["end"]
    canonical_orfs.loc[m_reverse_canonical, "orf_end"] = canonical_orfs.loc[
        m_reverse_canonical, "start"
    ]

    m_reverse_extracted = extracted_orfs["strand"] == "-"
    extracted_orfs["orf_end"] = extracted_orfs["end"]
    extracted_orfs.loc[m_reverse_extracted, "orf_end"] = extracted_orfs.loc[
        m_reverse_extracted, "start"
    ]

    # then, find extracted ORFs with the same "orf_end" (and seqname, strand) as canonical ORFs
    merge_fields = ["seqname", "strand", "orf_end"]
    canonical_extracted_orf_ends = canonical_orfs.merge(
        extracted_orfs, on=merge_fields, suffixes=["_canonical", "_extracted"]
    )

    # finally, pull this into a set
    zip_it = zip(
        canonical_extracted_orf_ends["id_canonical"],
        canonical_extracted_orf_ends["id_extracted"],
    )
    canonical_extracted_matching_ends = {(c, a) for c, a in zip_it}

    msg = "Finding ORFs which exactly overlap the canonical ORFs"
    logger.info(msg)

    exact_matches = bed_utils.get_bed_overlaps(
        canonical_orf_exons, extracted_orf_exons, min_a_overlap=1, min_b_overlap=1
    )

    exact_match_orf_ids = {m.b_info for m in exact_matches}

    m_exact_orf_matches = extracted_orf_exons["id"].isin(exact_match_orf_ids)
    extracted_orf_exons = extracted_orf_exons[~m_exact_orf_matches]

    m_canonical = extracted_orfs["id"].isin(exact_match_orf_ids)
    label = labels["canonical"]
    extracted_orfs.loc[m_canonical, "orf_type"] = label

    msg = f"Found {len(exact_match_orf_ids)} canonical ORFs labeled as {label}"
    logger.info(msg)

    msg = "Finding truncated canonical ORFs"
    logger.info(msg)

    truncated_matches = bed_utils.get_bed_overlaps(
        canonical_orf_exons, extracted_orf_exons, min_b_overlap=1
    )

    truncated_match_ids = {
        m.b_info
        for m in truncated_matches
        if (m.a_info, m.b_info) in canonical_extracted_matching_ends
    }

    m_truncated_matches = extracted_orf_exons["id"].isin(truncated_match_ids)
    extracted_orf_exons = extracted_orf_exons[~m_truncated_matches]

    m_canonical_truncated = extracted_orfs["id"].isin(truncated_match_ids)

    msg = "Finding extended canonical ORFs"
    logger.info(msg)

    extended_matches = bed_utils.get_bed_overlaps(
        canonical_orf_exons, extracted_orf_exons, min_a_overlap=1
    )

    # For standard assembly, we also need to make sure that
    # all extended matches are fully contained within the
    # transcript structure (i.e start upstream but otherwise
    # have the same structure).
    if args.nonoverlapping_label is None:

        transcript_matches = bed_utils.get_bed_overlaps(
            annotated_exons, extracted_orf_exons, min_b_overlap=1
        )
        transcript_match_pairs = {(m.a_info, m.b_info) for m in transcript_matches}

        extended_match_ids = {
            m.b_info
            for m in extended_matches
            if (m.a_info, m.b_info) in transcript_match_pairs
            and (m.a_info, m.b_info) in canonical_extracted_matching_ends
        }

    else:

        extended_match_ids = {
            m.b_info
            for m in extended_matches
            if (m.a_info, m.b_info) in canonical_extracted_matching_ends
        }

    m_extended_matches = extracted_orf_exons["id"].isin(extended_match_ids)
    extracted_orf_exons = extracted_orf_exons[~m_extended_matches]

    m_canonical_extended = extracted_orfs["id"].isin(extended_match_ids)
    m_canonical_variants = m_canonical_truncated | m_canonical_extended

    label = labels["{}canonical_variant".format(args.label_prefix)]
    extracted_orfs.loc[m_canonical_variants, "orf_type"] = label

    msg = (
        f"Found {len(extended_match_ids | truncated_match_ids)} "
        f"canonical_variant ORFs labeled as {label}"
    )
    logger.info(msg)

    msg = (
        "Finding within canonical ORFs that do not share an "
        "annotated stop codon with a canonical ORF (e.g. in "
        "frame stop, out-of-frame)"
    )
    logger.info(msg)

    within_ids = {
        m.b_info for m in truncated_matches if m.b_info not in truncated_match_ids
    }

    m_within_matches = extracted_orf_exons["id"].isin(within_ids)
    extracted_orf_exons = extracted_orf_exons[~m_within_matches]

    m_within = extracted_orfs["id"].isin(within_ids)
    label = labels["{}internal".format(args.label_prefix)]
    extracted_orfs.loc[m_within, "orf_type"] = label

    msg = f"Found {len(within_ids)} internal ORFs labeled as {label}"
    logger.info(msg)

    # find all overlapping ORFs
    msg = "Finding all UTR overlap matches"
    logger.info(msg)
    out_of_frame_matches = bed_utils.get_bed_overlaps(
        canonical_orf_exons, extracted_orf_exons
    )

    leader_matches = bed_utils.get_bed_overlaps(five_prime_exons, extracted_orf_exons)

    trailer_matches = bed_utils.get_bed_overlaps(three_prime_exons, extracted_orf_exons)

    msg = (
        "Labeling ORFs which have (out-of-frame) overlaps with both a "
        "canonical ORF and annotated leaders or trailers"
    )
    logger.info(msg)

    # We need to choose how to ensure that up-/downstream overlaps are unique.
    # Where an ORF overlaps both the 5'UTR and the 3'UTR of different same
    # sense overlapping transcripts, it is assigned by default to the downstream overlap.
    # For de novo, everything is labeled as overlap.

    leader_match_pairs = {(m.a_info, m.b_info) for m in leader_matches}
    trailer_match_pairs = {(m.a_info, m.b_info) for m in trailer_matches}

    if args.nonoverlapping_label is None:

        # For standard assembly, we also need to make sure that
        # all overlap matches are fully contained within the
        # transcript structure.
        transcript_matches = bed_utils.get_bed_overlaps(
            annotated_exons, extracted_orf_exons, min_b_overlap=1
        )

        transcript_match_pairs = {(m.a_info, m.b_info) for m in transcript_matches}

        leader_overlap_pairs = {
            (m.a_info, m.b_info)
            for m in out_of_frame_matches
            if (m.a_info, m.b_info) in leader_match_pairs
            and (m.a_info, m.b_info) not in trailer_match_pairs
            and (m.a_info, m.b_info) in transcript_match_pairs
        }

        trailer_overlap_pairs = {
            (m.a_info, m.b_info)
            for m in out_of_frame_matches
            if (m.a_info, m.b_info) in trailer_match_pairs
            and (m.a_info, m.b_info) not in leader_match_pairs
            and (m.a_info, m.b_info) in transcript_match_pairs
        }

        # We do not assign preference where the ORF overlaps both sides
        # of the coding sequence on the same transcript, any ORF
        # satisfying both will be labeled simply as overlap.
        overlap_ids = {
            m.b_info
            for m in out_of_frame_matches
            if (m.a_info, m.b_info) in leader_match_pairs
            and (m.a_info, m.b_info) in trailer_match_pairs
            and (m.a_info, m.b_info) in transcript_match_pairs
        }

        trailer_overlap_ids = {
            pair[1] for pair in trailer_overlap_pairs if pair[1] not in overlap_ids
        }

        leader_overlap_ids = {
            pair[1]
            for pair in leader_overlap_pairs
            if pair[1] not in trailer_overlap_ids and pair[1] not in overlap_ids
        }

        m_overlap_matches = extracted_orf_exons["id"].isin(overlap_ids)
        extracted_orf_exons = extracted_orf_exons[~m_overlap_matches]

        m_leader_overlap_matches = extracted_orf_exons["id"].isin(leader_overlap_ids)
        extracted_orf_exons = extracted_orf_exons[~m_leader_overlap_matches]

        m_five_prime_overlap = extracted_orfs["id"].isin(leader_overlap_ids)
        label = labels["{}five_prime_overlap".format(args.label_prefix)]
        extracted_orfs.loc[m_five_prime_overlap, "orf_type"] = label

        msg = (
            f"Found {len(leader_overlap_ids)} five_prime_overlap ORFs "
            f"labeled as {label}"
        )
        logger.info(msg)

        m_trailer_overlap_matches = extracted_orf_exons["id"].isin(trailer_overlap_ids)
        extracted_orf_exons = extracted_orf_exons[~m_trailer_overlap_matches]

        m_three_prime_overlap = extracted_orfs["id"].isin(trailer_overlap_ids)
        label = labels["{}three_prime_overlap".format(args.label_prefix)]
        extracted_orfs.loc[m_three_prime_overlap, "orf_type"] = label

        msg = (
            f"Found {len(trailer_overlap_ids)} three_prime_overlap ORFs "
            f"labeled as {label}"
        )
        logger.info(msg)

    else:

        overlap_ids = {m.b_info for m in out_of_frame_matches}
        overlap_ids |= {m.b_info for m in leader_matches}
        overlap_ids |= {m.b_info for m in trailer_matches}

        m_overlap_matches = extracted_orf_exons["id"].isin(overlap_ids)
        extracted_orf_exons = extracted_orf_exons[~m_overlap_matches]

    m_overlap = extracted_orfs["id"].isin(overlap_ids)
    label = labels["{}overlap".format(args.label_prefix)]
    extracted_orfs.loc[m_overlap, "orf_type"] = label

    msg = f"Found {len(overlap_ids)} overlap ORFs labeled as {label}"
    logger.info(msg)

    msg = "Finding ORFs completely within 5' or 3' leaders"
    logger.info(msg)

    leader_matches = bed_utils.get_bed_overlaps(
        five_prime_exons, extracted_orf_exons, min_b_overlap=1
    )

    leader_ids = {m.b_info for m in leader_matches}

    m_leader_matches = extracted_orf_exons["id"].isin(leader_ids)
    extracted_orf_exons = extracted_orf_exons[~m_leader_matches]

    m_five_prime = extracted_orfs["id"].isin(leader_ids)
    label = labels["{}five_prime".format(args.label_prefix)]
    extracted_orfs.loc[m_five_prime, "orf_type"] = label

    msg = f"Found {len(leader_ids)} five_prime ORFs labeled as {label}"
    logger.info(msg)

    trailer_matches = bed_utils.get_bed_overlaps(
        three_prime_exons, extracted_orf_exons, min_b_overlap=1
    )

    trailer_ids = {m.b_info for m in trailer_matches}

    m_trailer_matches = extracted_orf_exons["id"].isin(trailer_ids)
    extracted_orf_exons = extracted_orf_exons[~m_trailer_matches]

    m_three_prime = extracted_orfs["id"].isin(trailer_ids)
    label = labels["{}three_prime".format(args.label_prefix)]
    extracted_orfs.loc[m_three_prime, "orf_type"] = label

    msg = f"Found {len(trailer_ids)} three_prime ORFs labeled as {label}"
    logger.info(msg)

    msg = "Finding ORFs completely within annotated, non-coding transcripts"
    logger.info(msg)

    noncoding_matches = bed_utils.get_bed_overlaps(
        noncoding_exons, extracted_orf_exons, min_b_overlap=1
    )

    noncoding_ids = {m.b_info for m in noncoding_matches}

    m_noncoding_matches = extracted_orf_exons["id"].isin(noncoding_ids)
    extracted_orf_exons = extracted_orf_exons[~m_noncoding_matches]

    m_noncoding = extracted_orfs["id"].isin(noncoding_ids)
    label = labels["{}noncoding".format(args.label_prefix)]
    extracted_orfs.loc[m_noncoding, "orf_type"] = label

    msg = f"Found {len(noncoding_ids)} noncoding ORFs labeled as {label}"
    logger.info(msg)

    # all of the remaining ORFs fall into the "suspect" category
    suspect_ids = {orf_id for orf_id in extracted_orf_exons["id"]}

    m_suspect = extracted_orfs["id"].isin(suspect_ids)
    label = labels["{}suspect".format(args.label_prefix)]
    extracted_orfs.loc[m_suspect, "orf_type"] = label

    n_suspect_ids = len(suspect_ids)
    msg = f"Remaining {n_suspect_ids} ORFs labeled as {label}"
    logger.info(msg)

    m_no_orf_type = extracted_orfs["orf_type"].isnull()
    msg = "Found {} unlabeled ORFs".format(sum(m_no_orf_type))
    logger.info(msg)

    msg = "Writing ORFs and labels to disk"
    logger.info(msg)

    extracted_orfs = bed_utils.sort(extracted_orfs)

    additional_columns = ["orf_num", "orf_len"]
    fields = bed_utils.bed12_field_names + additional_columns
    orfs_genomic = extracted_orfs[fields]
    bed_utils.write_bed(orfs_genomic, args.extracted_orfs)

    label_columns = ["orf_num", "id", "orf_type", "transcripts"]
    extracted_orfs = extracted_orfs[label_columns]
    bed_utils.write_bed(extracted_orfs, args.out)


if __name__ == "__main__":
    main()
