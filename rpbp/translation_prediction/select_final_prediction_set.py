#! /usr/bin/env python3

import argparse
import logging
import pandas as pd

import Bio.Seq

import bio_utils.bio as bio
import bio_utils.bed_utils as bed_utils
import bio_utils.fastx_utils as fastx_utils
import misc.logging_utils as logging_utils
import misc.utils as utils
import misc.parallel as parallel

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

def get_best_overlapping_orf(merged_ids, predicted_orfs):
    if len(merged_ids) == 1:
        m_id = predicted_orfs['id'] == merged_ids[0]
        return predicted_orfs[m_id].iloc[0]
    
    m_orfs = predicted_orfs['id'].isin(merged_ids)
    
    sorted_orfs = predicted_orfs[m_orfs].sort_values('bayes_factor_mean', ascending=False)
    best_overlapping_orf_id = sorted_orfs.iloc[0]
    return best_overlapping_orf_id

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

    parser.add_argument('--select-longest-by-stop', help="If this flag is given, then "
        "the selected ORFs will be merged based on stop codons. In particular, only the "
        "longest translated ORF at each stop codon will be selected.", action='store_true')

    parser.add_argument('--select-best-overlapping', help="If this flag is given, then "
        "only the ORF with the highest estimated Bayes factor will be kept among each "
        "set of overlapping ORFs. N.B. This filter is applied *AFTER* selecting the "
        "longest ORF at each stop codon, if the --select-longest-by-stop flag is "
        "given.", action='store_true')

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
        
    
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    # first, extract all of the predictions which exceed the threshold
    msg = "Reading Bayes factor information"
    logger.info(msg)
    
    bayes_factors = bed_utils.read_bed(args.bayes_factors)

    if args.filter_non_canonical_overlaps:
        args.filtered_orf_types.extend(non_canonical_overlap_orf_types)

    if len(args.filtered_orf_types) > 0:
        filtered_orf_types_str = ','.join(args.filtered_orf_types)
        msg = "Filtering these ORF types: {}".format(filtered_orf_types_str)
        logger.info(msg)

        m_orf_types = bayes_factors['orf_type'].isin(args.filtered_orf_types)
        bayes_factors = bayes_factors[~m_orf_types]

    msg = "Identifying ORFs which meet the prediction thresholds"
    logger.info(msg)

    all_orfs, bf_orfs, chisq_orfs = ribo_utils.get_predicted_orfs(
        bayes_factors, 
        min_bf_mean=args.min_bf_mean, 
        max_bf_var=args.max_bf_var, 
        min_bf_likelihood=args.min_bf_likelihood,
        min_length=args.min_length,
        chisq_alpha=args.chisq_significance_level,
        select_longest_by_stop=args.select_longest_by_stop
    )
    
    if args.use_chi_square:
        predicted_orfs = chisq_orfs
    else:
        predicted_orfs = bf_orfs

    msg = "Number of selected ORFs: {}".format(len(predicted_orfs))
    logger.info(msg)

    if args.select_best_overlapping:

        msg = "Finding overlapping ORFs"
        logger.info(msg)

        merged_intervals = bed_utils.merge_all_intervals(predicted_orfs)

        msg = "Selecting best among overlapping ORFs"
        logger.info(msg)

        predicted_orfs = parallel.apply_iter_simple(
            merged_intervals['merged_ids'], 
            get_best_overlapping_orf, 
            predicted_orfs, 
            progress_bar=True
        )

        predicted_orfs = pd.DataFrame(predicted_orfs)

    msg = "Sorting selected ORFs"
    logger.info(msg)

    predicted_orfs = bed_utils.sort(predicted_orfs)

    msg = "Writing selected ORFs to disk"
    logger.info(msg)
    bed_utils.write_bed(predicted_orfs, args.predicted_orfs)

    # now get the sequences
    msg = "Extracting predicted ORFs DNA sequence"
    logger.info(msg)

    split_exons = True
    transcript_sequences = bed_utils.get_all_bed_sequences(
        predicted_orfs, 
        args.fasta, 
        split_exons
    )

    fastx_utils.write_fasta(transcript_sequences, args.predicted_dna_sequences,
        compress=False)

    # translate the remaining ORFs into protein sequences
    msg = "Converting predicted ORF sequences to amino acids"
    logger.info(msg)

    records = fastx_utils.get_read_iterator(args.predicted_dna_sequences)
    protein_records = {
        r[0]: Bio.Seq.translate(r[1]) for r in records
    }
    
    fastx_utils.write_fasta(
        protein_records.items(), 
        args.predicted_protein_sequences, 
        compress=False
    )

if __name__ == '__main__':
    main()

