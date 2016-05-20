#! /usr/bin/env python3

import argparse
import logging

import misc.bio as bio
import misc.utils as utils

default_types = []

canonical_types = ['canonical_truncated', 'canonical', 'canonical_extended']

non_coding_types = ['five_prime', 'three_prime', 'noncoding', 'novel']

non_canonical_types = ['five_prime_overlap', 'suspect_overlap', 'five_prime', 'three_prime', 
    'three_prime_overlap', 'noncoding', 'novel_suspect_overlap', 'novel']

extended_canonical_types = ['canonical_truncated', 'canonical', 'canonical_extended', 
    'five_prime_overlap', 'suspect_overlap', 'three_prime_overlap',
    'novel_suspect_overlap']

all_types = ['canonical', 'within', 'canonical_extended', 'canonical_truncated',
       'three_prime_overlap', 'suspect_overlap', 'noncoding', 'five_prime',
       'five_prime_overlap', 'three_prime', 'novel_suspect_overlap',
       'novel']

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script extracts ORFs of specified types from any of the BED12+ "
        "files created by this pipeline. The 'type' arguments will be aggregated for "
        "the final filter.")
    parser.add_argument('bed', help="The input BED12+ file")
    parser.add_argument('out', help="The filtered BED12+ file")
    parser.add_argument('--standard', help="If this flag is present, then only the "
        "standard BED12 fieds will be included in the output file", action='store_true')
    parser.add_argument('--compress', help="if this flag is present, then the "
        "output will be gzipped (but the output filename will not be changed",
        action='store_true')
    parser.add_argument('--types', help="The types can be explicitly listed if the "
        "predefined types are not appropriate", nargs='*', default=default_types)
    parser.add_argument('--canonical', help="Give this flag to include the canonical "
        "(plus extended and truncated) types", action='store_true')
    parser.add_argument('--noncoding', help="Give this flag to include the noncoding "
        "types: 'five_prime', 'three_prime', 'noncoding', 'novel'", action='store_true')
    parser.add_argument('--non-canonical', help="Give this flag to include all types "
        "except the canonical (extended and truncated) types", action='store_true')
    parser.add_argument('--extended-canonical', help="Give this flag to include all "
        "types which overlap canonical coding transcripts in some way", action='store_true')
    parser.add_argument('--all', help="Give this flag to select all types. This is "
        "mostly helpful to fix the BED score and color columns", action='store_true')

    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    types = args.types

    if args.canonical:
        types = types + canonical_types

    if args.noncoding:
        types = types + non_coding_types

    if args.non_canonical:
        types = types + non_canonical_types

    if args.extended_canonical:
        types = types + extended_canonical_types

    if args.all:
        types = all_types

    types_str = '\n'.join(types)
    msg = "Extracting the following types:\n{}".format(types_str)
    logging.info(msg)

    bed_df = bio.read_bed(args.bed)
    bed_df['score'] = 0
    bed_df['color'] = 0
    m_types = bed_df['orf_type'].isin(types)

    fields = bed_df.columns
    if args.standard:
        fields = bio.bed12_field_names

    bio.write_bed(bed_df.loc[m_types, fields], args.out, compress=args.compress)

if __name__ == '__main__':
    main()
