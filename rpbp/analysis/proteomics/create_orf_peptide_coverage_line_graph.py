#! /usr/bin/env python3

import matplotlib
matplotlib.rc('text', usetex=True)

import argparse
import logging
import re

import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd

import misc.parallel as parallel
import misc.utils as utils

default_num_procs = 1

default_min_length = 0

default_title = ""
default_fontsize = 20
default_note_fontsize = 15
default_line_width = 4

def get_orf_coverage(orf):
    """ This function calculates the percentage of an ORF covered by the
        matching peptides.
        
        It is adopted from a StackOverflow answer:
        http://stackoverflow.com/questions/4173904/string-coverage-optimization-in-python
    
    """
    
    peptides = orf['peptide_matches']
    peptides = peptides.replace(";", "|")
    pat = re.compile('(' + peptides + ')')
    
    orf_sequence = orf['orf_sequence']
    orf_id = orf['orf_id']
    num_matches = orf['num_matches']
    peptide_length = len(orf_sequence)
    
    sout = list(orf_sequence)
    i = 0
    match = pat.search(orf_sequence)
    while match:
        span = match.span()
        sout[span[0]:span[1]] = ['x']*(span[1]-span[0])
        i = span[0]+1
        match = pat.search(orf_sequence, i)
    covered_orf_sequence = ''.join(sout)
    
    num_covered_positions = covered_orf_sequence.count('x')
    coverage = num_covered_positions / peptide_length
    
    ret = {
        'orf_id': orf_id,
        'num_matches': num_matches,
        'covered_orf_sequence': covered_orf_sequence,
        'num_covered_positions': num_covered_positions,
        'coverage': coverage,
        'peptide_length': peptide_length
    }
    
    return pd.Series(ret)


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script creates a plot showing the fraction of predicted ORFs "
        "which have a set amount of peptide coverage.")
    parser.add_argument('rpbp_peptide_matches', help="The (csv) file containing the peptides "
        "matching to each ORF predicted as translated using Rp-Bp (produced by "
        "get-orf-peptide-matches)")
    parser.add_argument('rpchi_peptide_matches', help="The (csv) file containing the peptides "
        "matching to each ORF predicted as translated using Rp-chi (produced by "
        "get-orf-peptide-matches)")
    parser.add_argument('out', help="The output (image) file")

    parser.add_argument('-l', '--min-length', help="The minimum length for ORFs (in "
        "nucleotides) to consider in the analyis", type=int, default=default_min_length)

    parser.add_argument('-p', '--num-procs', help="The number of processors to use", type=int,
        default=default_num_procs)

    parser.add_argument('--title', default=default_title)
    parser.add_argument('--fontsize', type=int, default=default_fontsize)
    parser.add_argument('--note-fontsize', type=int, default=default_note_fontsize)
    parser.add_argument('--line-width', type=int, default=default_line_width)

    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    msg = "Reading predictions"
    logging.info(msg)
    rpbp_peptide_matches = pd.read_csv(args.rpbp_peptide_matches)
    rpchi_peptide_matches = pd.read_csv(args.rpchi_peptide_matches)

    if args.min_length > 0:
        msg = "Filtering predictions by: length > {}".format(args.min_length)
        logging.warning(msg)

        # multiply by 3 because the orf sequences are amino acid sequences
        bf_lengths = rpbp_peptide_matches['orf_sequence'].str.len() * 3
        m_bf_length = bf_lengths > args.min_length
        rpbp_peptide_matches = rpbp_peptide_matches[m_bf_length]


        chisq_lengths = rpchi_peptide_matches['orf_sequence'].str.len() * 3
        m_chisq_length = chisq_lengths > args.min_length
        rpchi_peptide_matches = rpchi_peptide_matches[m_chisq_length]
        
    msg = "Calculating Rp-Bp coverage"
    logging.info(msg)
    bf_coverage = parallel.apply_parallel(rpbp_peptide_matches, 
                                            args.num_procs, 
                                            get_orf_coverage, 
                                            progress_bar=True)
    bf_coverage = pd.DataFrame(bf_coverage)

    msg = "Calculating Rp-chi coverage"
    logging.info(msg)
    chisq_coverage = parallel.apply_parallel(rpchi_peptide_matches, 
                                                args.num_procs,
                                                get_orf_coverage, 
                                                progress_bar=True)
    chisq_coverage = pd.DataFrame(chisq_coverage)

    msg = "Creating image"
    logging.info(msg)

    # plot the empirical distribution of ORF lengths
    hist_min = 0
    hist_max = 1.1
    hist_step = 0.05
    hist_range = (hist_min, hist_max)
    hist_bins = np.arange(hist_min, hist_max, hist_step)

    bf_covered_hist, b = np.histogram(bf_coverage['coverage'], 
                                        bins=hist_bins, 
                                        range=hist_range, 
                                        density=True)
    
    chisq_covered_hist, b = np.histogram(chisq_coverage['coverage'], 
                                            bins=hist_bins, 
                                            range=hist_range, 
                                            density=True)

    # now, normalize the histograms
    bf_covered_hist = bf_covered_hist / np.sum(bf_covered_hist)
    chisq_covered_hist = chisq_covered_hist / np.sum(chisq_covered_hist)

    # multiply by 100 to give actual percentages
    bf_covered_hist = 100 * bf_covered_hist
    chisq_covered_hist = 100 * chisq_covered_hist

    hist_bins = 100*hist_bins

    fig, ax = plt.subplots(figsize=(10,5))

    cm = plt.cm.gist_earth

    x = np.arange(len(bf_covered_hist))

    bf_label = r'\textsc{Rp-Bp}'
    ax.plot(x, 
            bf_covered_hist, 
            color=cm(0.1), 
            label=bf_label, 
            linewidth=args.line_width, 
            linestyle='--', 
            marker='^')

    chisq_label = r'\textsc{Rp-$\chi^2$}'
    ax.plot(x, 
            chisq_covered_hist, 
            color=cm(0.3), 
            label=chisq_label, 
            linewidth=args.line_width,
            linestyle='-.',
            marker='D')

    ax.set_xlabel('Peptide Coverage (\%)', fontsize=args.fontsize)
    ax.set_ylabel('\% of predicted ORFs', fontsize=args.fontsize)
    
    if args.title is not None and len(args.title) > 0:
        ax.set_title(args.title, fontsize=args.fontsize)

    # only show every 20% on the x-axis
    ax.set_xticks(x[::4])
    ax.set_xticklabels(hist_bins[::4])

    def my_formatter_fun(x, p):
        return "${:d}$".format(20*p)
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(my_formatter_fun))

    # hide the "0" tick label
    yticks = ax.yaxis.get_major_ticks()
    yticks[0].label1.set_visible(False)

    ax.set_xlim((0, len(bf_covered_hist)-1))
    ax.set_ylim((0,10))

    ax.legend(loc='upper right', fontsize=args.fontsize)
    ax.tick_params(axis='both', which='major', labelsize=args.note_fontsize)

    fig.savefig(args.out, bbox_inches='tight')

if __name__ == '__main__':
    main()
