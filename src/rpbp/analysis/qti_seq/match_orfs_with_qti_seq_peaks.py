#! /usr/bin/env python3

import argparse
import pandas as pd
import logging

import pybedtools

import pbiotools.utils.bio as bio
import pbiotools.misc.utils as utils

logger = logging.getLogger(__name__)

default_seqname_prefix = ""
default_output_prefix = "qti_"


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script matches QTI-seq peaks from (Gao et al., 2015) to ORFs "
        "based on genomic coordinates.",
    )
    parser.add_argument("orfs", help="The ORFs (bed12+) file")
    parser.add_argument("qti_peaks", help="The QTI-seq peak (BED6) files")
    parser.add_argument("out", help="The augmented ORFs (BED12+) file")

    parser.add_argument(
        "--output-prefix",
        help="A string to prefix before all of the "
        "fields related to the closest QTI-seq peak (if there is one)",
        default=default_output_prefix,
    )

    parser.add_argument(
        "--seqname-prefix",
        help="If present, this string is prepended "
        "to all of the ORF seqnames. It is then removed again in the final output.",
        default=default_seqname_prefix,
    )

    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    programs = ["closestBed"]
    utils.check_programs_exist(programs)

    msg = "Reading ORFs"
    logger.info(msg)

    orfs = bio.read_bed(args.orfs)

    # we need to keep a copy that we use later for output
    orfs_copy = orfs.copy()

    # for matching qti-seq peaks, we only want to consider the start position of each ORF

    # for forward strand ORFs, replace orf_genomic_end with orf_genomic_start
    msg = "Updating genomic positions to consider only start codon"
    logger.info(msg)
    m_forward = orfs["strand"] == "+"
    # forward_orfs = orfs[mask_forward]
    # forward_orfs['end'] = forward_orfs['start'] + 1
    orfs.loc[m_forward, "end"] = orfs.loc[m_forward, "start"] + 1

    # for reverse ORFs, replace orf_genomic_start with orf_genomic_end
    m_reverse = orfs["strand"] == "-"
    # reverse_orfs = orfs[mask_reverse]
    # reverse_orfs['start'] = reverse_orfs['end'] - 1
    orfs.loc[m_reverse, "start"] = orfs.loc[m_reverse, "end"] - 1

    # join together the orf start positions, correct the seqname and sort for bedtools
    msg = "Converting ORF data frame to pybedtools"
    logger.info(msg)
    # orfs_start_only = pd.concat([forward_orfs, reverse_orfs])
    # orfs_start_only['seqname'] = args.seqname_prefix + orfs_start_only['seqname']
    # orfs_start_only = orfs_start_only.sort_values(['seqname', 'start'])
    # orfs_bed = pybedtools.BedTool.from_dataframe(orfs_start_only)

    orfs["seqname"] = args.seqname_prefix + orfs["seqname"]
    orfs = orfs.sort_values(["seqname", "start"])
    orfs_bed = pybedtools.BedTool.from_dataframe(orfs)

    msg = "Reading QTI peaks"
    qti_bed_df = bio.read_bed(args.qti_peaks)
    qti_bed_df.columns = [
        "{}_{}".format(args.output_prefix, c) for c in qti_bed_df.columns
    ]

    # and covert to bed
    msg = "Converting QTI peaks data frame to pybedtools"
    logger.info(msg)

    chr_field = "{}_chr".format(args.output_prefix)
    start_field = "{}_start".format(args.output_prefix)
    qti_bed_df = qti_bed_df.sort_values([chr_field, start_field])
    qti_bed = pybedtools.BedTool.from_dataframe(qti_bed_df)

    msg = "Finding closest QTI peak for all ORFs"
    logger.info(msg)
    # s means to consider strandedness
    # D means to report the distance
    # a means to report upstream positions (relative to orfs_start) as negative
    closest_bed = orfs_bed.closest(qti_bed, s=True, D="a")

    # covert back to a df for clean up
    msg = "Converting closest results back to data frame"
    logger.info(msg)

    peak_distance_field = "{}_peak_distance".format(args.output_prefix)
    closest_bed_fields = (
        list(orfs.columns) + list(qti_bed_df.columns) + [peak_distance_field]
    )
    closest_df = closest_bed.to_dataframe(names=closest_bed_fields, index_col=False)

    # now join the relevant fields back to the original ORF data frame
    fields_to_join = list(qti_bed_df.columns) + [peak_distance_field, "id"]
    closest_df = closest_df[fields_to_join]

    msg = "Joining closest results to original ORFs"
    logger.info(msg)
    orf_qti_df = pd.merge(orfs_copy, closest_df, on="id", how="left")
    orf_qti_df = orf_qti_df.sort_values(["seqname", "start"])

    # and write this out as a bed12+ file
    msg = "Writing joined BED12+ file to disk"
    logger.info(msg)

    bio.write_bed(orf_qti_df, args.out)


if __name__ == "__main__":
    main()
