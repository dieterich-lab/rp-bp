#! /usr/bin/env python3

import argparse

import collections
import pandas as pd

import bio_utils.bed_utils as bed_utils
import bio_utils.mygene_utils as mygene_utils
import misc.logging_utils as logging_utils
import misc.parallel as parallel
import misc.utils as utils

import logging
import misc.logging_utils as logging_utils
logger = logging.getLogger(__name__)

default_scopes = []

fields_to_keep = [
       'orf_len',
       'bayes_factor_mean',
       'bayes_factor_var',
       'x_1_sum',
       'x_2_sum',
       'x_3_sum'
]

fields_to_keep = bed_utils.bed12_field_names + fields_to_keep

orf_id_info = collections.namedtuple(
    "orf_id_info",
    "transcript_id,seqname,strand,start,end,length,orf_id"
)
def parse_orf_id(orf_id, trim=True):
    transcript_id, s = orf_id.split("_")
    seqname, s, strand = s.split(":")
    start, end, length = bed_utils.parse_exon_start_end_length(s)
    
    if trim:
        transcript_id = transcript_id.split(".")[0]
    
    ret = orf_id_info(
        transcript_id,
        seqname,
        strand,
        start,
        end,
        length,
        orf_id
    )
    
    return ret

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script uses the mygene.info service to find annotations "
        "for the transcripts associated with the ORFs in the given bed file. In "
        "particular, it extracts information from Swiss-Prot, TrEMBL, Interpro, "
        "PDB, Pfam, PROSITE, the Gene Ontology, and KEGG.")

    parser.add_argument('bed', help="The bed file")
    parser.add_argument('out', help="The output file. Its type will be inferred "
        "from its extension.")

    parser.add_argument('--do-not-trim', help="By default, the script will "
        "attempt to trim transcript identifiers such that they are valid Ensembl "
        "identifiers. If this flag is given, no trimming will take place.",
        action='store_true')

    parser.add_argument('--scopes', help="A list of scopes to use when querying "
        "mygene.info. Please see the documentation for more information about "
        "valid scopes: http://mygene.info/doc/query_service.html#available_fields",
        nargs='*', default=default_scopes)

    parser.add_argument('--do-not-convert-ids', help="By default, the script will "
        "treat the identifiers in the file as transcript identifiers. It first "
        "maps those to gene identifiers, and then it uses those to find the "
        "gene annotations. If the identifiers are already gene ids (or whatever "
        "is specified by scopes), then the first mapping is not necessary and "
        "can be skipped using this flag.", action='store_true')
    
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    convert_ids = not args.do_not_convert_ids

    msg = "Reading the bed file"
    logger.info(msg)
    bed = bed_utils.read_bed(args.bed)
    bed = bed[fields_to_keep]

    msg = "Extracting transcript ids"
    logger.info(msg)
    trim = not args.do_not_trim
    orf_ids = parallel.apply_iter_simple(bed['id'], parse_orf_id, trim)
    orf_ids_df = pd.DataFrame(orf_ids)

    if convert_ids:
        msg = "Querying transcript to gene id mapping"
        logger.info(msg)
        gene_ids = mygene_utils.get_transcript_to_gene_mapping(orf_ids_df['transcript_id'])
    else:
        gene_ids = pd.DataFrame()
        gene_ids['transcript_id'] = orf_ids_df['transcript_id']
        gene_ids['gene_id'] = orf_ids_df['transcript_id']

    msg = "Querying gene annotations"
    logger.info(msg)
    res_df = mygene_utils.query_mygene(gene_ids['gene_id'])

    msg = "Combining gene annotations with transcript ids"
    logger.info(msg)
    res_df = gene_ids.merge(res_df, on='gene_id', how='inner')

    msg = "Combining transcript annotations with ORF ids"
    logger.info(msg)
    orf_ids_fields = ['transcript_id', 'orf_id']
    res_df = orf_ids_df[orf_ids_fields].merge(res_df, on='transcript_id', how='inner')

    msg = "Combining ORF annotations with ORF predictions"
    logger.info(msg)
    res_df = bed.merge(res_df, left_on='id', right_on='orf_id', how='left')

    msg = "Writing ORF annotations to disk"
    logger.info(msg)
    utils.write_df(res_df, args.out, index=False)

if __name__ == '__main__':
    main()
