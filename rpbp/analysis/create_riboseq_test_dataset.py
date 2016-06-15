#! /usr/bin/env python3

import argparse
import itertools
import logging
import yaml

import pysam

import misc.bio as bio
import misc.utils as utils
import rpbp.filenames as filenames

default_reference = 'I'
default_max_reads = 100000

def get_first_token(s):
    return s.split()[0]

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script extracts reads of various types from a processed dataset "
        "to create an \"interesting\" test dataset.\n\nN.B. This script is not "
        "particularly efficient and is not intended for external use.")
    parser.add_argument('config', help="The (yaml) config file")
    parser.add_argument('name', help="The name of the dataset to use to create the test data")
    parser.add_argument('out', help="The output (fasta.gz) which contains reads of various "
        "types, subject to the other parameters")

    parser.add_argument('-r', '--reference', help="The name of the reference sequence (chromosome) "
        "from which aligned reads will be extracted", default=default_reference)

    parser.add_argument('-m', '--max-reads', help="At most <max_reads> reads of each type "
        "will be included in the final output", type=int, default=default_max_reads)
    
    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    config = yaml.load(open(args.config))

    note = config.get('note', None)

    ###
    msg = "Reading alignments from BAM file"
    logging.info(msg)

    bam_file = filenames.get_riboseq_bam(config['riboseq_data'], 
        args.name, is_unique=True, note=note)
    bam = pysam.AlignmentFile(bam_file)
    
    alignments = bam.fetch(reference=args.reference)
    num_alignments = bam.count(reference=args.reference)
    alignment_qnames = {get_first_token(a.qname) for a in alignments}

    ###
    msg = "Extracting a similar number of rRNA reads"
    logging.info(msg)

    with_rrna = filenames.get_with_rrna_fastq(config['riboseq_data'], args.name, note=note)
        
    rrna = bio.get_read_iterator(with_rrna, is_fasta=False)
    rrna = itertools.islice(rrna, num_alignments)
    rrna_qnames = {get_first_token(read[0]) for read in rrna}

    ###
    msg = "Extracting a similar number of reads which do not uniquely map to the genome"
    logging.info(msg)

    # first, pull out the qnames of all alignments
    all_alignments = bam.fetch()
    all_alignment_qnames = {get_first_token(a.qname) for a in all_alignments}

    # iterate over all reads which passed the rRNA and quality filtering
    without_rrna_file = filenames.get_without_rrna_fastq(config['riboseq_data'], args.name, note=note)
    without_rrna = bio.get_read_iterator(without_rrna_file, is_fasta=False)
    without_rrna_qnames = {get_first_token(read[0]) for read in without_rrna}

    no_mapping_qnames = without_rrna_qnames - all_alignment_qnames

    ###
    msg = "Extracting a similar number of reads which are filtered due to quality issues"
    logging.info(msg)

    # first, pull in all the reads and their names

    msg = "Reading all reads into a dictionary"
    logging.debug(msg)

    raw_data_file = config['riboseq_samples'][args.name]
    raw_data = bio.get_fasta_dict(raw_data_file, is_fasta=False, key_fn=get_first_token)
    raw_data_qnames = set(raw_data.keys())

    msg = "Reading quality scores into dictionary"
    logging.debug(msg)

    raw_data_qual = bio.get_fastq_qual_dict(raw_data_file, key_fn=get_first_token)

    # now, the reads which _did_ pass quality filtering
    msg = "Reading reads which pass quality filtering into a set"
    logging.debug(msg)

    without_adapters_file = filenames.get_without_adapters_fastq(config['riboseq_data'], args.name, note=note)
    without_adapters = bio.get_read_iterator(without_adapters_file, is_fasta=False)
    without_adapters_qnames = {get_first_token(read[0]) for read in without_adapters}

    # and pull out the qnames of the reads which did not pass quality filtering
    filtered_reads_qnames = raw_data_qnames - without_adapters_qnames

    ###
    msg = "Constructing the set of reads to output"
    logging.info(msg)
        
    alignment_raw_data = {qname: raw_data[qname] 
        for qname in itertools.islice(alignment_qnames, args.max_reads)}

    rrna_raw_data = {qname: raw_data[qname] 
        for qname in itertools.islice(rrna_qnames, args.max_reads)}

    no_mapping_raw_data = {qname: raw_data[qname] 
        for qname in itertools.islice(no_mapping_qnames, args.max_reads)}

    filtered_reads_raw_data = {qname: raw_data[qname] 
        for qname in itertools.islice(filtered_reads_qnames, args.max_reads)}

    out_raw_data = alignment_raw_data
    out_raw_data.update(rrna_raw_data)
    out_raw_data.update(no_mapping_raw_data)
    out_raw_data.update(filtered_reads_raw_data)

    ###
    msg = "Writing sequences to disk"
    logging.info(msg)

    bio.write_fastq(out_raw_data, raw_data_qual, args.out, progress_bar=True)

if __name__ == '__main__':
    main()
