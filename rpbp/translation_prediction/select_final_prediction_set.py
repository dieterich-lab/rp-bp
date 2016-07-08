#! /usr/bin/env python3

import argparse
import logging
import pandas as pd

import pybedtools
from Bio.Seq import Seq

import misc.bio as bio
import misc.utils as utils

import riboutils.ribo_utils as ribo_utils

logger = logging.getLogger(__name__)

default_min_bf_mean = 5
default_min_bf_likelihood = 0.5
default_max_bf_var = None
default_min_length = 20

default_chisq_significance_level = 0.01

default_chi_square_field = 'chi_square_p'
default_bf_mean_field = 'bayes_factor_mean'
default_bf_var_field = 'bayes_factor_var'
default_profile_sum_field = 'profile_sum'

default_filtered_orf_types = []
non_canonical_overlap_orf_types = [
    'five_prime_overlap',
    'suspect_overlap',
    'three_prime_overlap',
    'within'
]
non_canonical_overlap_orf_types_str = ','.join(non_canonical_overlap_orf_types)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Given a list of ORFs with associated Bayes factors and a fasta "
        "sequence file, this script extracts the sequences of the ORFs whose Bayes factor "
        "exceeds the given threshold. Finally, biopython is used to translate the "
        "selected ORFs into protein sequences.\n\n"
        "The min-length and minimum-profile-sum filters are applied in the obvious way.\n\n"
        "For both BF and chi-square predictions, only ORFs which have more reads in the "
        "first reading frame than either of the other two will be selected as translated. "
        "(This is called the 'frame filter' below.)\n\n"
        "The selection based on Bayes factors follows this logic: if max_bf_var is given, "
        "then it and min_bf_mean are taken as a hard threshold on the estimated Bayes "
        "factor mean. If min_bf_likelihood is given, then this min_bf_mean is taken as the "
        "boundary value; that is, an ORF is \"translated\" if:\n\n"
        "\t\t[P(bf > min_bf_mean)] > min_bf_likelihood\n\n"
        "If both max_bf_var and min_bf_likelihood are None, then min_bf_mean is taken as a "
        "hard threshold on the mean for selecting translated ORFs.\n\n"
        "If both max_bf_var and min_bf_likelihood are given, then both filters will be "
        "applied and the result will be the intersection.\n\n"
        "If the --use-chi-square option is given, the significance value is "
        "Bonferroni-corrected based on the number of ORFs which meet the length, profile "
        "and frame filters.")

    parser.add_argument('bayes_factors', help="The file containing the ORFs and Bayes' "
        "factors (BED12+)")
    parser.add_argument('fasta', help="The *genome* fasta file")
    parser.add_argument('predicted_orfs', help="The (output) BED12+ file containing "
        "the predicted ORFs.")
    parser.add_argument('predicted_dna_sequences', help="The (output) fasta file "
        "containing the predicted ORF sequences, as DNA sequences")
    parser.add_argument('predicted_protein_sequences', help="The (output) fasta file "
        "containing the predicted ORF sequences, as protein sequences")

    parser.add_argument('--min-length', help="The minimum length to predict an ORF "
        "as translated", type=int, default=default_min_length)
    
    parser.add_argument('--min-bf-mean', help="The minimum Bayes' factor mean to predict "
        "an ORF as translated (use --help for more details)", 
        type=float, default=default_min_bf_mean)
    parser.add_argument('--max-bf-var', help="The maximum Bayes' factor variance to predict "
        "an ORF as translated (use --help for more details)", 
        type=float, default=default_max_bf_var)

    parser.add_argument('--min-bf-likelihood', help="If given, then this is taken a threshold "
        "on the likelihood of translation (use --help for more details)", 
        type=float, default=default_min_bf_likelihood)
   
    parser.add_argument('--use-chi-square', help="If this flag is present, the the "
        "chi square value will be used to predict ORFs rather than the Bayes' factor",
        action='store_true')
    parser.add_argument('--chisq-significance-level', help="If using chi square, then this "
        "value is Bonferroni corrected and used as the significance cutoff", type=float, 
        default=default_chisq_significance_level)

    parser.add_argument('--filtered-orf-types', help="A list of ORF types which will be "
        "removed before selecting the final prediction set.", nargs='*', 
        default=default_filtered_orf_types)

    parser.add_argument('--filter-non-canonical-overlaps', help="If this flag is given, then "
        "--filtered-orf-types will be extended with the non-canonical overlap types ({})."
        .format(non_canonical_overlap_orf_types_str), action='store_true')
        
    
    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    programs = ['fastaFromBed']
    utils.check_programs_exist(programs)

    # first, extract all of the predictions which exceed the threshold
    msg = "Reading Bayes factor information"
    logger.info(msg)
    
    bayes_factors = bio.read_bed(args.bayes_factors)

    if args.filter_non_canonical_overlaps:
        args.filtered_orf_types.extend(non_canonical_overlap_orf_types)

    if len(args.filtered_orf_types) > 0:
        filtered_orf_types_str = ','.join(args.filtered_orf_types)
        msg = "Filtering these ORF types: {}".format(filtered_orf_types_str)
        logger.info(msg)

        m_filtered_orf_types = bayes_factors['orf_type'].isin(args.filtered_orf_types)
        bayes_factors = bayes_factors[~m_filtered_orf_types]

    msg = "Identifying ORFs which meet the prediction thresholds"
    longest_orfs, bf_orfs, chisq_orfs = ribo_utils.get_predicted_orfs(bayes_factors, 
                                                       min_bf_mean=args.min_bf_mean, 
                                                       max_bf_var=args.max_bf_var, 
                                                       min_bf_likelihood=args.min_bf_likelihood,
                                                       min_length=args.min_length,
                                                       chisq_alpha=args.chisq_significance_level)
    
    if args.use_chi_square:
        longest_predicted_orfs = chisq_orfs
    else:
        longest_predicted_orfs = bf_orfs

    msg = "Number of longest ORFs: {}".format(len(longest_predicted_orfs))
    logger.info(msg)

    bio.write_bed(longest_predicted_orfs, args.predicted_orfs)

    # now get the sequences
    msg = "Extracting predicted ORFs DNA sequence"
    logger.info(msg)

    # first, convert the dataframe to a bedtool
    filtered_predictions_bed = pybedtools.BedTool.from_dataframe(longest_predicted_orfs[bio.bed12_field_names])

    # now, pull out the sequences
    predicted_dna_sequences = filtered_predictions_bed.sequence(fi=args.fasta, 
            split=True, s=True, name=True)
    predicted_dna_sequences.save_seqs(args.predicted_dna_sequences)

    # translate the remaining ORFs into protein sequences
    msg = "Converting predicted ORF sequences to amino acids"
    logger.info(msg)

    records = bio.get_read_iterator(args.predicted_dna_sequences)
    protein_records = {
        r[0]: str(Seq(r[1]).translate()) for r in records
    }
    bio.write_fasta(protein_records, args.predicted_protein_sequences, compress=False)

if __name__ == '__main__':
    main()

