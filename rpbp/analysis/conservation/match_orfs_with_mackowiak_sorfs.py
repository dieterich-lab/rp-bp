#! /usr/bin/env python3

import argparse
import logging
import pandas as pd

import pybedtools

import misc.bio as bio
import misc.utils as utils

logger = logging.getLogger(__name__)

default_seqname_prefix = ""
default_output_prefix = "mackowiak"

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script matches sORFs from (Mackowiak et al., 2015) to ORFs "
        "based on genomic coordinates.")
    parser.add_argument('orfs', help="The ORFs (BED12+) file")
    parser.add_argument('mackowiak_orfs', help="The ORFs given BED12 format. This should "
        "be the output from the mackowiak2bed.py script (or that output subsequently "
        "processed); it should *not* be the raw supplementary files.")
    parser.add_argument('out', help="The labeled ORFs (bed.gz) file")
    
    parser.add_argument('--seqname-prefix', help="If present, this string is prepended "
        "to all of the ORF seqnames. It is then removed again in the final output.", 
        default=default_seqname_prefix)

    parser.add_argument('--output-prefix', help="A prefix to add to the fields related to "
        "the matching Mackowiak ORFs in the output", default=default_output_prefix)

    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    programs = ['intersectBed']
    utils.check_programs_exist(programs)

    msg = "Reading Mackowiak ORFs and converting to pybed"
    logger.info(msg)

    mackowiak_sorfs = bio.read_bed(args.mackowiak_orfs)
    mackowiak_sorfs = mackowiak_sorfs[bio.bed12_field_names]
    mackowiak_sorfs.columns = ["{}_{}".format(args.output_prefix, c)
                                    for c in mackowiak_sorfs.columns]
    mackowiak_sorfs_bed = pybedtools.BedTool.from_dataframe(mackowiak_sorfs)

    msg = "Reading our ORFs and converting to pybed"
    logger.info(msg)

    orfs = bio.read_bed(args.orfs)

    if len(args.seqname_prefix) > 0:
        orfs['seqname'] = args.seqname_prefix + orfs['seqname']

    orfs_bed = pybedtools.BedTool.from_dataframe(orfs)

    # now, find our ORFs which intersect any of mackowiak's orfs
    msg = "Intersecting ORFs and converting to data frame"
    logger.info(msg)

    # u means to report "a" features which intersect any "b" feature
    # wa means to write the "a" bed features
    # wb means to write the (intersecting) "b" bed features
    # split means to consider the "blocks"
    # loj means "left outer join", with "a" (orfs_bed) considered as left
    # s means to consider strandedness
    loj_intersection = orfs_bed.intersect(mackowiak_sorfs_bed, split=True, wa=True, wb=True, loj=True, s=True)

    # convert the intersection results back to a data frame
    intersection_bed_fields = list(orfs.columns) + list(mackowiak_sorfs.columns)
    loj_df = loj_intersection.to_dataframe(names=intersection_bed_fields)

    if len(args.seqname_prefix) > 0:
        pattern = "^{}".format(args.seqname_prefix)
        loj_df['seqname'] = loj_df['seqname'].str.replace(pattern, "")

    # and write to disk
    msg = "Writing intersected ORFs to disk as BED12+ file"
    logger.info(msg)

    bio.write_bed(loj_df, args.out)

if __name__ == '__main__':
    main()
