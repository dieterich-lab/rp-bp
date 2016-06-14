#! /usr/bin/env python3

import argparse
import yaml
import logging
import os
import pandas as pd

import misc.latex as latex
import misc.utils as utils
import rpbp.filenames as filenames

default_tmp = None

default_project_name = "ribosome profiling data"
default_min_metagene_profile_count = 1000
default_min_metagene_profile_bayes_factor_mean = 5
default_max_metagene_profile_bayes_factor_var = 5
default_min_visualization_count = 500
default_num_procs = 1

abstract = """
This document shows the results of preprocessing steps. 
In particular, it shows the amount of reads filtered at each step in the pipeline. 
Additionally, it includes ``metagene'' profile figures for all of the samples.
"""

read_filtering_label = "fig:mapping-info"
length_distribution_section_label = "sec:read-length-distribution"
periodicity_label = "sec:periodicity"

mapping_and_filtering_text = """
Our \\riboseq processing pipeline consists of the following steps.

\\begin{enumerate}
\\item Adapters and low quality reads are removed.
\\item Reads mapping to rRNA are removed.
\\item Reads are aligned to the genome.
\\item Reads with multiple alignments are removed.
\\item Reads with length that do not result in a strong periodic signal are removed.
\\item Remaing reads are used to construct the filtered genome profile
\\end{enumerate}

Figure~\\ref{fig:mapping-info}(left) shows the number of reads remaining after each stage in our preprocessing pipeline for all samples.
Figure~\ref{fig:mapping-info}(right) shows a ``zoomed in'' version which does not include the reads of poor quality and that mapped to ribosomal sequences.
"""

read_length_distribution_text = """
This section shows the distribution of read lengths \\textbf{for reads which uniquely map to the genome}.
"""

read_filtering_caption = """
The number of reads lost at each stage in the mapping pipeline. 
\texttt{wrong\_length} refers to reads which have a length that does not result in a strong periodic signal (see Section~\ref{sec:periodicity}).
"""

def create_figures(config_file, config, name, offsets_df, args):
    """ This function creates all of the figures in the preprocessing report
        for the given dataset.
    """
    note = config.get('note', None)
    # visualize the metagene profiles
    msg = "{}: Visualizing metagene profiles and Bayes' factors".format(name)
    logging.info(msg)
    logging.debug("overwrite: {}".format(args.overwrite))

    metagene_profiles = filenames.get_metagene_profiles(config['riboseq_data'], 
        name, is_unique=True, note=note)
    
    profile_bayes_factor = filenames.get_metagene_profiles_bayes_factors(config['riboseq_data'],
        name, is_unique=True, note=note)

    mp_df = pd.read_csv(metagene_profiles)

    min_read_length = int(offsets_df['length'].min())
    max_read_length = int(offsets_df['length'].max())

    for length in range(min_read_length, max_read_length+1):

        mask_length = offsets_df['length'] == length

        # TODO: it is not clear why, but it seems sometimes no rows match
        if sum(mask_length) == 0:
            continue
        length_row = offsets_df[mask_length].iloc[0]
               
        # make sure we have enough reads to visualize
        if length_row['highest_peak_profile_sum'] < args.min_visualization_count:
            continue
        
        # visualize the metagene profile
        title = "Periodicity, {}, length {}".format(name, length)
        metagene_profile_image = filenames.get_riboseq_metagene_profile_image(config['riboseq_data'], 
            name, image_type='eps', is_unique=True, length=length, note=note)
        cmd = ("visualize-metagene-profile {} {} {} --title \"{}\"".format(
            metagene_profiles, length, metagene_profile_image, title))
        in_files = [metagene_profiles]
        out_files = [metagene_profile_image]
        utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=True)

        # and the Bayes' factor
        title = "Metagene profile Bayes' factors, {}, length {}".format(name, length)
        metagene_profile_image = filenames.get_metagene_profile_bayes_factor_image(
            config['riboseq_data'], name, image_type='eps', is_unique=True, length=length, note=note)
        cmd = ("visualize-metagene-profile-bayes-factor {} {} {} --title \"{}\" "
            "--font-size 25".format(profile_bayes_factor, length, metagene_profile_image, title))
        in_files = [profile_bayes_factor]
        out_files = [metagene_profile_image]
        utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=True)

def create_read_filtering_plots(config_file, config, args):
        
    # get the filtering counts
    note = config.get('note', None)
    read_filtering_counts = filenames.get_riboseq_read_filtering_counts(config['riboseq_data'], note=note)
    overwrite_str = ""
    if args.overwrite:
        overwrite_str = "--overwrite"

    tmp_str = ""
    if args.tmp is not None:
        tmp_str = "--tmp {}".format(args.tmp)

    logging_str = utils.get_logging_options_string(args)

    procs_str = "--num-procs {}".format(args.num_procs)
    cmd = "get-all-read-filtering-counts {} {} {} {} {} {}".format(config_file, 
        read_filtering_counts, overwrite_str, procs_str, tmp_str, logging_str)
    in_files = [config_file]
    out_files = [read_filtering_counts]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=True)

    # and visualize them
    read_filtering_image = filenames.get_riboseq_read_filtering_counts_image(config['riboseq_data'], note=note)
    title = "Read filtering counts"
    cmd = "visualize-read-filtering-counts {} {} --title \"{}\"".format(read_filtering_counts, 
        read_filtering_image, title)
    in_files = [read_filtering_counts]
    out_files=[read_filtering_image]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=True)

    # and visualize the filtering without the rrna
    n = "no-rrna-{}".format(note)
    read_filtering_image = filenames.get_riboseq_read_filtering_counts_image(config['riboseq_data'], note=n)
    title = "Read filtering counts, no ribosomal matches"
    cmd = "visualize-read-filtering-counts {} {} --title \"{}\" --without-rrna".format(
        read_filtering_counts, read_filtering_image, title)

    in_files = [read_filtering_counts]
    out_files=[read_filtering_image]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=True)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script creates a simple latex document containing the read "
        "filtering images, metagene profiles and analysis, and standard section text.")
    parser.add_argument('config', help="The (yaml) config file for the project")
    parser.add_argument('out', help="The path for the output files")

    parser.add_argument('--num-procs', help="The number of processors to use for counting "
        "the mapped reads. This is only useful up to the number of samples in the project.",
        type=int, default=default_num_procs)
    parser.add_argument('--overwrite', help="If this flag is present, existing files will "
        "be overwritten.", action='store_true')

    parser.add_argument('--min-visualization-count', help="Read lengths with fewer than this "
        "number of reads will not be included in the report.", type=int, 
        default=default_min_visualization_count)

    parser.add_argument('--tmp', help="Intermediate files (such as fastqc reports when "
        "they are first generated) will be written here", default=default_tmp)

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

    note = config.get('note', None)


    # make sure the path to the output file exists
    os.makedirs(args.out, exist_ok=True)

   
    # first, create the read filtering information
    create_read_filtering_plots(args.config, config, args)

    min_metagene_profile_count = config.get(
        "min_metagene_profile_count", default_min_metagene_profile_count)

    min_metagene_profile_bayes_factor_mean = config.get(
        "min_metagene_profile_bayes_factor_mean", default_min_metagene_profile_bayes_factor_mean)

    max_metagene_profile_bayes_factor_var = config.get(
        "max_metagene_profile_bayes_factor_var", default_max_metagene_profile_bayes_factor_var)


    project_name = config.get("project_name", default_project_name)
    title = "Preprocessing results for {}".format(project_name)
   
    header = latex.get_header_text(title, abstract)
    footer = latex.get_footer_text()

    tex_file = os.path.join(args.out, "preprocessing-report.tex")
    with open(tex_file, 'w') as out:

        out.write(header)
        out.write("\n")

        latex.section(out, "Introduction")

        latex.clearpage(out)
        latex.newpage(out)

        latex.section(out, "Mapping and filtering")
        out.write(mapping_and_filtering_text)

        # the read filtering figures

        read_filtering_image = filenames.get_riboseq_read_filtering_counts_image(
            config['riboseq_data'], note=note)
    
        n = "no-rrna-{}".format(note)
        no_rrna_read_filtering_image = filenames.get_riboseq_read_filtering_counts_image(
            config['riboseq_data'], note=n)

        latex.begin_figure(out)
        latex.write_graphics(out, read_filtering_image, width=0.75)
        latex.write_graphics(out, no_rrna_read_filtering_image, width=0.75)
        latex.write_caption(out, read_filtering_caption, label=read_filtering_label)
        latex.end_figure(out)

        # the read length distributions
        latex.section(out, "Read length distributions", label=length_distribution_section_label)
        out.write(read_length_distribution_text)

        msg = "Writing length distribution figures"
        logging.info(msg)

        for name, data in config['riboseq_samples'].items():
            unique_filename_fastqc_read_lengths = filenames.get_riboseq_bam_fastqc_read_lengths(
                config['riboseq_data'], name, is_unique=True, note=note)

            caption = "Read length distribution of mapped reads for {}".format(name)
            latex.begin_figure(out)
            latex.write_graphics(out, unique_filename_fastqc_read_lengths, width=0.5)
            latex.write_caption(out, caption)
            latex.end_figure(out)


        latex.section(out, "Periodicity", label=periodicity_label)

        for name, data in config['riboseq_samples'].items():
            msg = "Processing sample: {}".format(name)
            logging.info(msg)

            logging.debug("overwrite: {}".format(args.overwrite))
    
            periodic_offsets = filenames.get_periodic_offsets(config['riboseq_data'], 
                name, is_unique=True, note=note)
            offsets_df = pd.read_csv(periodic_offsets)

            min_read_length = int(offsets_df['length'].min())
            max_read_length = int(offsets_df['length'].max())
    
            create_figures(args.config, config, name, offsets_df, args)


            # first, the read length figure
            without_rrna_fastqc_read_lengths = filenames.get_riboseq_bam_fastqc_read_lengths(
                config['riboseq_data'], name, note=note)
            read_lengths_caption = "Read length distribution of mapped read for {}".format(name)

            latex.begin_figure(out)
            latex.write_graphics(out, without_rrna_fastqc_read_lengths, width=0.5)
            latex.write_caption(out, read_lengths_caption)
            latex.end_figure(out)

            latex.begin_figure(out)

            i = 0
            for length in range(min_read_length, max_read_length + 1):
                msg = "Processing length: {}".format(length)
                logging.info(msg)

                # check which offset is used
                
                # select the row for this length
                mask_length = offsets_df['length'] == length

                # TODO: this is sometimes length 0. why?
                if sum(mask_length) == 0:
                    continue

                length_row = offsets_df[mask_length].iloc[0]

                # now, check all of the filters
                offset = int(length_row['highest_peak_offset'])
                offset_status = "Used for analysis"
                
                if length_row['highest_peak_bf_mean'] < min_metagene_profile_bayes_factor_mean:
                    offset_status = "BF mean too small"

                if length_row['highest_peak_bf_var'] > max_metagene_profile_bayes_factor_var:
                    offset_status = "BF variance too high"

                if length_row['highest_peak_profile_sum'] < min_metagene_profile_count:
                    offset_status = "Count too small"
                
                if length_row['highest_peak_profile_sum'] < args.min_visualization_count:
                    msg = "Not enough reads of this length. Skipping."
                    logging.warning(msg)
                    continue

                out.write(name.replace('_', '-'))
                out.write(", length: ")
                out.write(str(length))
                out.write(", offset: ")

                out.write(str(offset))

                out.write(", status: ")
                out.write(offset_status)

                out.write("\n\n")

                metagene_profile_image = filenames.get_riboseq_metagene_profile_image(
                    config['riboseq_data'], name, image_type='eps', is_unique=True, length=length, note=note)
                
                bayes_factor_image = filenames.get_metagene_profile_bayes_factor_image(
                    config['riboseq_data'], name, image_type='eps', is_unique=True, length=length, note=note)

                latex.write_graphics(out, metagene_profile_image, width=0.6)
                latex.write_graphics(out, bayes_factor_image, width=0.39)
                                
                i += 1

                if i % 5 == 0:
                    latex.write_caption(out, name)
                    latex.end_figure(out)
                    latex.clearpage(out)
                    latex.begin_figure(out)

            latex.write_caption(out, name)
            latex.end_figure(out)

        out.write(footer)

    os.chdir(args.out)
    cmd = "pdflatex -shell-escape preprocessing-report"
    utils.check_call(cmd)
    utils.check_call(cmd) # call again to fix references

if __name__ == '__main__':
    main()

