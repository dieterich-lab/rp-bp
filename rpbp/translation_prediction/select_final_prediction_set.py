#! /usr/bin/env python3

import argparse
import logging
import pandas as pd

import pybedtools
from Bio.Seq import Seq

import misc.bio as bio
import misc.utils as utils

import rpbp.rpbp_utils as rpbp_utils

default_min_bf_mean = 5
default_max_bf_var = 2
default_min_length = 20
default_minimum_profile_sum = 5

default_chisq_significance_level = 0.01

default_chi_square_field = 'chi_square_p'
default_bf_mean_field = 'bayes_factor_mean'
default_bf_var_field = 'bayes_factor_var'
default_profile_sum_field = 'profile_sum'

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Given a list of ORFs with associated Bayes' factors and a fasta "
        "sequence file, this script extracts the sequences of the ORFs whose Bayes' factor "
        "exceeds the given threshold. Finally, biopython is used to translate the "
        "selected ORFs into protein sequences.")
    parser.add_argument('bayes_factors', help="The file containing the ORFs and Bayes' "
        "factors (BED12+)")
    parser.add_argument('fasta', help="The *genome* fasta file")
    parser.add_argument('predicted_orfs', help="The (output) BED12+ file containing "
        "the predicted ORFs.")
    parser.add_argument('predicted_dna_sequences', help="The (output) fasta file "
        "containing the predicted ORF sequences, as DNA sequences")
    parser.add_argument('predicted_protein_sequences', help="The (output) fasta file "
        "containing the predicted ORF sequences, as protein sequences")

    parser.add_argument('--all', help="If this flag is present, then all ORF sequences "
        "will be extracted.", action='store_true')


    parser.add_argument('--min-bf-mean', help="The minimum Bayes' factor mean to predict "
        "an ORF as translated", type=float, default=default_min_bf_mean)
    parser.add_argument('--max-bf-var', help="The maximum Bayes' factor variance to predict "
        "an ORF as translated", type=float, default=default_max_bf_var)
    parser.add_argument('--min-length', help="The minimum length to predict an ORF "
        "as translated", type=int, default=default_min_length)
    parser.add_argument('--minimum-profile-sum', help="The minimum read count (or however "
        "the profile sum is calculated) to consider for prediction", type=float,
        default=default_minimum_profile_sum)

    parser.add_argument('--use-chi-square', help="If this flag is present, the the "
        "chi square value will be used to predict ORFs rather than the Bayes' factor",
        action='store_true')
    parser.add_argument('--chisq-alpha', help="If using chi square, then this "
        "value is Bonferroni corrected and used as the significance cutoff", type=float, 
        default=default_chisq_significance_level)
    
    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    programs = ['fastaFromBed']
    utils.check_programs_exist(programs)

    # first, extract all of the predictions which exceed the threshold
    msg = "Identifying ORFs which meet the prediction thresholds"
    logging.info(msg)

    bayes_factors = bio.read_bed(args.bayes_factors)

    longest_orfs, bf_orfs, chisq_orfs = rpbp_utils.get_predicted_orfs(bayes_factors, 
                                                       min_signal=args.minimum_profile_sum, 
                                                       min_bf_mean=args.min_bf_mean, 
                                                       max_bf_var=args.max_bf_var, 
                                                       min_length=args.min_length,
                                                       chisq_alpha=args.chisq_alpha)
    
    if args.use_chi_square:
        longest_predicted_orfs = chisq_orfs
    else:
        longest_predicted_orfs = bf_orfs

    msg = "Number of longest ORFs: {}".format(len(longest_predicted_orfs))
    logging.info(msg)

    bio.write_bed(longest_predicted_orfs, args.predicted_orfs)

    # now get the sequences
    msg = "Extracting predicted ORFs DNA sequence"
    logging.info(msg)

    # first, convert the dataframe to a bedtool
    filtered_predictions_bed = pybedtools.BedTool.from_dataframe(longest_predicted_orfs[bio.bed12_field_names])

    # now, pull out the sequences
    predicted_dna_sequences = filtered_predictions_bed.sequence(fi=args.fasta, 
            split=True, s=True, name=True)
    predicted_dna_sequences.save_seqs(args.predicted_dna_sequences)

    # translate the remaining ORFs into protein sequences
    msg = "Converting predicted ORF sequences to amino acids"
    logging.info(msg)

    records = bio.get_read_iterator(args.predicted_dna_sequences)
    protein_records = {
        r[0]: str(Seq(r[1]).translate()) for r in records
    }
    bio.write_fasta(protein_records, args.predicted_protein_sequences, compress=False)

if __name__ == '__main__':
    main()

