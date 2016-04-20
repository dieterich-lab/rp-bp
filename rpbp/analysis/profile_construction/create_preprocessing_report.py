#! /usr/bin/env python3

import argparse
import yaml
import logging
import os
import pandas as pd

import misc.utils as utils
import riboseq_utils.riboseq_filenames as filenames

def create_figures(config_file, config, name, overwrite):
    """ This function creates all of the figures in the preprocessing report
        for the given dataset.
    """
    # visualize the metagene profiles
    msg = "{}: Visualizing metagene profiles and Bayes' factors".format(name)
    logging.info(msg)

    logging.debug("overwrite: {}".format(overwrite))

    seqids_to_keep_str = utils.get_config_argument(config, 'seqids_to_keep')
    read_length_range = range(config['min_read_length'], config['max_read_length']+1)

    periodicity_bayes_factor = filenames.get_riboseq_metagene_periodicity_bayes_factors(config['riboseq_data'],
        name, is_unique=True)

    for length in read_length_range:
        metagene_profile = filenames.get_riboseq_metagene_profile(config['riboseq_data'], 
            name, is_unique=True, length=length)
        
        # visualize the metagene profile
        title = "Periodicity, {}, length {}".format(name, length)
        series_label = "Length {}".format(length)
        metagene_profile_image = filenames.get_riboseq_metagene_profile_image(config['riboseq_data'], 
            name, image_type='eps', is_unique=True, length=length)
        cmd = ("visualize-metagene-profile {} {} --title '{}' --series-label '{}'".format(
            metagene_profile, metagene_profile_image, title, series_label))
        in_files = [metagene_profile]
        out_files = [metagene_profile_image]
        utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=overwrite, call=True)

        # and the Bayes' factor
        title = "Metagene profile Bayes' factors, {}, length {}".format(name, length)
        series_label = "Length {}".format(length)
        metagene_profile_image = filenames.get_riboseq_metagene_profile_bayes_factor_image(config['riboseq_data'], 
            name, image_type='eps', is_unique=True, length=length)
        cmd = ("visualize-metagene-profile-bayes-factor {} {} {} --title \"{}\" --series-label \"{}\" "
            "--font-size 25".format(periodicity_bayes_factor, length, metagene_profile_image, 
            title, series_label))
        in_files = [metagene_profile]
        out_files = [metagene_profile_image]
        utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=overwrite, call=True)

 
        
    # get the filtering counts
    msg = "{}: Getting reading filtering counts".format(name)
    logging.info(msg)

    read_filtering_counts = filenames.get_riboseq_read_filtering_counts(config['riboseq_data'])
    cmd = "get-all-read-filtering-counts {} {}".format(config_file, read_filtering_counts)
    in_files = [config_file]
    out_files = [read_filtering_counts]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=overwrite, call=True)

    # and visualize them
    
    read_filtering_counts = filenames.get_riboseq_read_filtering_counts(config['riboseq_data'])
    overwrite_str = ""
    if overwrite:
        overwrite_str = "--overwrite"
    cmd = "get-all-read-filtering-counts {} {} {}".format(config_file, read_filtering_counts, overwrite_str)
    in_files = [config_file]
    out_files = [read_filtering_counts]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=overwrite, call=True)

    read_filtering_image = filenames.get_riboseq_read_filtering_counts_image(config['riboseq_data'])
    title = "Read filtering counts"
    cmd = "visualize-read-filtering-counts {} {} --title \"{}\"".format(read_filtering_counts, read_filtering_image, title)
    in_files = [read_filtering_counts]
    out_files=[read_filtering_image]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=overwrite, call=True)

    # and visualize the filtering without the rrna
    read_filtering_image = filenames.get_riboseq_read_filtering_counts_image(config['riboseq_data'], note="no-rrna")
    title = "Read filtering counts, no ribosomal matches"
    cmd = "visualize-read-filtering-counts {} {} --title \"{}\" --without-rrna".format(read_filtering_counts, read_filtering_image, title)

    in_files = [read_filtering_counts]
    out_files=[read_filtering_image]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=overwrite, call=True)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script creates a simple latex document containing the read "
        "filtering images, metagene profiles and analysis, and standard section text.")
    parser.add_argument('config', help="The (yaml) config file for the project")
    parser.add_argument('out', help="The latex output file")
    parser.add_argument('--overwrite', help="If this flag is present, existing files will "
        "be overwritten.", action='store_true')

    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    config = yaml.load(open(args.config))

    programs =  [   'visualize-metagene-profile',
                    'visualize-metagene-profile-bayes-factor',
                    'get-all-read-filtering-counts',
                    'fastqc',
                    'java',
                    'samtools',
                    'visualize-read-filtering-counts'
                ]
    utils.check_programs_exist(programs)


    # make sure the path to the output file exists
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    header, footer = get_header_and_footer_text()


    with open(args.out, 'w') as out:

        out.write(header)
        out.write("\n")

        for name in config['sample_names']:
            msg = "Processing sample: {}".format(name)
            logging.info(msg)

            logging.debug("overwrite: {}".format(args.overwrite))
            create_figures(args.config, config, name, args.overwrite)

            out.write("\\begin{figure}\n")
            out.write("\t\\centering\n")

            best_periodicity_and_offsets = filenames.get_riboseq_best_periodicity_and_offsets(config['riboseq_data'], name, is_unique=True)
            best_periodicity_and_offsets = pd.read_csv(best_periodicity_and_offsets)
            
            i = 0
            for length in range(config['min_read_length'], config['max_read_length'] + 1):
                msg = "Processing length: {}".format(length)
                logging.info(length)

                out.write(name.replace('_', '-'))
                out.write(", length: ")
                out.write(str(length))
                out.write(", selected offset: ")

                # check which offset is used
                
                # select the row for this length
                mask_length = best_periodicity_and_offsets['length'] == length
                length_row = best_periodicity_and_offsets[mask_length].iloc[0]

                offset = "BF too small"
                if length_row['largest_count_bf'] > config['min_periodicity_bf']:
                    offset = int(length_row['largest_count_offset'])

                if length_row['largest_count'] < config['min_periodicity_count']:
                    offset = "Count too small"
                
                out.write(str(offset))

                out.write("\n\n")

                metagene_profile_image = filenames.get_riboseq_metagene_profile_image(config['riboseq_data'], name, image_type='eps', is_unique=True, length=length)
                metagene_profile_image = metagene_profile_image.replace(".eps", "}.eps")
                
                bayes_factor_image = filenames.get_riboseq_metagene_profile_bayes_factor_image(config['riboseq_data'], name, image_type='eps', is_unique=True, length=length)
                bayes_factor_image = bayes_factor_image.replace(".eps", "}.eps")
                
                out.write("\t\\includegraphics[width=0.6\\textwidth,keepaspectratio]{{")
                out.write(metagene_profile_image)
                out.write("}\n")
                
                out.write("\t\\includegraphics[width=0.39\\textwidth,keepaspectratio]{{")
                out.write(bayes_factor_image)
                out.write("}\n")

                
                i += 1

                if i % 5 == 0:
                    out.write("\t\\caption{{{0}}}\n".format(name.replace("_", "-")))
                    out.write("\\end{figure}\n")
                    out.write("\n")
                    out.write("\\begin{figure}\n")
                    out.write("\t\\centering\n")

            out.write("\t\\caption{{{0}}}\n".format(name.replace("_", "-")))
            out.write("\\end{figure}\n")
            out.write("\n")
            

        out.write(footer)

if __name__ == '__main__':
    main()

def get_header_and_footer_text():
    header = r"""
\documentclass[a4paper,10pt]{article} % For LaTeX2e
\usepackage[utf8x]{inputenc}
\usepackage{amsmath}
\usepackage{graphicx}
\usepackage{framed}
\usepackage{float}

\usepackage{morefloats}

\usepackage{xspace}

\pagestyle{empty}
\setlength{\parskip}{2mm}
\setlength{\parindent}{0cm}




\title{Preprocessing results for Human S-tau 1 data}

\def\|#1{\ensuremath{\mathtt{#1}}}
\def\!#1{\ensuremath{\mathbf{#1}}}
\def\*#1{\ensuremath{\mathcal{#1}}}

\newcommand*{\defeq}{\mathrel{\vcenter{\baselineskip0.5ex \lineskiplimit0pt
                     \hbox{\scriptsize.}\hbox{\scriptsize.}}}%
                     =}


\DeclareMathOperator*{\argmax}{arg\,max}

\newtheorem{definition}{Definition}
\newtheorem{theorem}{Theorem}

\newcommand{\bitseq}{\textsc{BitSeq}\xspace}
\newcommand{\flexbar}{\textsc{Flexbar}\xspace}
\newcommand{\bowtiep}{\textsc{Bowtie2}\xspace}
\newcommand{\starp}{\textsc{STAR}\xspace}

\newcommand{\te}{\textsc{TE}\xspace}
\newcommand{\riboseq}{ribo-seq\xspace}
\newcommand{\rnaseq}{\textsc{RNA}-Seq\xspace}
\newcommand{\tets}{\ensuremath{\te_t^s}\xspace}
\newcommand{\ribots}{\ensuremath{ribo_t^s}\xspace}
\newcommand{\rnats}{\ensuremath{rna_t^s}\xspace}

\newcommand{\dep}{\ensuremath{\mathcal{M}^{dep}}}
\newcommand{\indep}{\ensuremath{\mathcal{M}^{indep}}}

\begin{document}


\maketitle

\begin{abstract}
This document shows the results of preprocessing steps. 
In particular, it shows the amount of reads filtered at each step in the pipeline. 
Additionally, it includes ``metagene'' periodicity figures for all of the samples.
\end{abstract}

\section{Introduction}


\section{Mapping and filtering}

Our \riboseq processing pipeline consists of the following steps.

\textbf{1}. Adapters and low quality reads are removed.\\
\textbf{2}. Reads mapping to rRNA are removed.\\
\textbf{3}. Reads are aligned to the genome.\\
\textbf{4}. Reads with multiple alignments are removed.\\
\textbf{5}. Reads with length that do not result in a strong periodic signal are removed.\\
\textbf{6}. Remaing reads are used to construct transcript signals.\\
\textbf{7}. Micropeptides, etc., are predicted using transcript signals.

Figure~\ref{fig:mapping-info}(left) shows the number of reads remaining after each stage in our preprocessing pipeline for all samples.
Figure~\ref{fig:mapping-info}(right) shows a ``zoomed in'' version which does not include the reads of poor quality and that mapped to ribosomal sequences.

\begin{figure}
    \centering
    \includegraphics[width=0.48\textwidth,keepaspectratio]{read-filtering-counts.eps}
    \includegraphics[width=0.48\textwidth,keepaspectratio]{{read-filtering-counts.no-rrna}.eps}
    \caption{The number of reads lost at each stage in the mapping pipeline. \texttt{wrong\_length} refers to reads which have a length that does not result in a strong periodic signal (see Section~\ref{sec:periodicity}). Please note the different colors in the figures. \label{fig:mapping-info}}
\end{figure}


\section{Periodicity}
\label{sec:periodicity}

        """

    footer = r"""

%\small
%\bibliographystyle{abbrv}
%\bibliography{library}

\end{document}
        """

    return (header, footer)
