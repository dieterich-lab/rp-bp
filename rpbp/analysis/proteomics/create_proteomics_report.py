#! /usr/bin/env python3

import argparse
import logging
import os
import yaml

import pbio.misc.latex as latex
import pbio.misc.utils as utils

import pbio.ribo.ribo_utils
import pbio.ribo.ribo_filenames as filenames


default_tmp = None
default_min_orf_length = 0
default_project_name = "riboseq proteomics report"
default_note = None
default_num_cpus = 1


def create_figures(config_file, config, name, args):
    """ This function creates all of the figures in the proteomics report
        for the given dataset.
    """

    out_note_str = config.get('note', None)
    if args.note is not None and len(args.note) > 0:
        out_note_str = args.note

    # keep multimappers?
    is_unique = not ('keep_riboseq_multimappers' in config)

    msg = "{}: Creating the ORF-peptide coverage line graph".format(name)
    logging.info(msg)

    try:
        lengths, offsets = riboutils.ribo_utils.get_periodic_lengths_and_offsets(
            config, name, is_unique=is_unique)
    except FileNotFoundError:
        msg = ("Could not parse out lengths and offsets for sample: {}. "
            "Skipping".format(name))
        logging.error(msg)
        return
    
    peptide_matches = filenames.get_riboseq_peptide_matches(
        config['riboseq_data'], name, length=lengths, offset=offsets, 
        is_unique=is_unique, note=out_note_str)
    
    chisq_peptide_matches = filenames.get_riboseq_peptide_matches(
        config['riboseq_data'], name, length=lengths, offset=offsets, 
        is_unique=is_unique, note=out_note_str, is_chisq=True)

    peptide_coverage_line_graph = filenames.get_peptide_coverage_line_graph(config['riboseq_data'], 
        name, length=lengths, offset=offsets, is_unique=is_unique, note=out_note_str)

    title_str = "--title {}".format(name)
    min_length_str = "--min-length {}".format(args.min_orf_length)
    num_cpus_str = "--num-cpus {}".format(args.num_cpus)

    cmd = "create-orf-peptide-coverage-line-graph {} {} {} {} {} {}".format(
        peptide_matches, chisq_peptide_matches, peptide_coverage_line_graph,
        title_str, min_length_str, num_cpus_str)

    in_files = [peptide_matches, chisq_peptide_matches]
    out_files = [peptide_coverage_line_graph]
    utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=True)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script creates a simple latex document.")
    parser.add_argument('config', help="The (yaml) config file for the project")
    parser.add_argument('out', help="The path for the output files")

    parser.add_argument('-l', '--min-orf-length', help="The minimum length for ORFs (in "
        "nucleotides) to consider in the analyis", type=int, default=default_min_orf_length)


    parser.add_argument('--num-cpus', help="The number of processors to use for counting "
        "the matching peptides coverage for each predicted ORF.",
        type=int, default=default_num_cpus)
    parser.add_argument('--overwrite', help="If this flag is present, existing files will "
        "be overwritten.", action='store_true')
        
    parser.add_argument('--note', help="If this option is given, it will be used in the "
        "filenames.\n\nN.B. This REPLACES the note in the config file.", default=default_note)


    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    config = yaml.load(open(args.config), Loader=yaml.FullLoader)

    # keep multimappers?
    is_unique = not ('keep_riboseq_multimappers' in config)

    programs =  [   
        'create-orf-peptide-coverage-line-graph'
    ]
    utils.check_programs_exist(programs)
    
    required_keys = [
        'riboseq_data',
        'riboseq_samples'
    ]
    utils.check_keys_exist(config, required_keys)

    # make sure the path to the output file exists
    os.makedirs(args.out, exist_ok=True)

    note_str = config.get('note', None)
    out_note_str = note_str

    if args.note is not None and len(args.note) > 0:
        out_note_str = args.note

    project_name = config.get("project_name", default_project_name)
    title = "Proteomics analysis for {}".format(project_name)
    abstract = "This document shows the results of proteomics analysis."

    header = latex.get_header_text(title, abstract)
    footer = latex.get_footer_text()

    tex_file = os.path.join(args.out, "proteomics-report.tex")

    with open(tex_file, 'w') as out:

        out.write(header)
        out.write("\n")

        title = "ORF peptide coverage"
        latex.section(out, title)


        for name, data in config['riboseq_samples'].items():
            msg = "Processing sample: {}".format(name)
            logging.info(msg)

            logging.debug("overwrite: {}".format(args.overwrite))

            create_figures(args.config, config, name, args)

            title = name
            latex.subsection(out, title)

            try:
                lengths, offsets = riboutils.ribo_utils.get_periodic_lengths_and_offsets(
                    config, name, is_unique=is_unique)
            except FileNotFoundError:
                msg = "Could not parse out lengths and offsets for sample: {}. Skipping".format(name)
                logging.error(msg)
                return
            
            peptide_coverage_line_graph = filenames.get_peptide_coverage_line_graph(config['riboseq_data'], 
                name, length=lengths, offset=offsets, is_unique=is_unique, note=out_note_str)

            if os.path.exists(peptide_coverage_line_graph):
                latex.begin_figure(out)
                latex.write_graphics(out, peptide_coverage_line_graph, width=0.9)
                latex.end_figure(out)

            else:
                text = "Problem creating ORF peptide coverage plot"
                out.write(text)
                out.write("\n")

            latex.clearpage(out)

        out.write(footer)

    
    os.chdir(args.out)
    cmd = "pdflatex -shell-escape proteomics-report"
    utils.check_call(cmd)
    utils.check_call(cmd) # call again to fix references            


if __name__ == '__main__':
    main()
