#! /usr/bin/env python3

import argparse
import yaml
import logging
import os
import pandas as pd
import shlex
import sys

import pbio.misc.latex as latex
import pbio.misc.logging_utils as logging_utils
import pbio.misc.parallel as parallel
import pbio.misc.shell_utils as shell_utils
import pbio.misc.slurm as slurm

import misc.utils as utils
import pbio.ribo.ribo_filenames as filenames
import pbio.ribo.ribo_utils as ribo_utils

logger = logging.getLogger(__name__)

default_note = None
default_tmp = None
default_image_type = 'eps'

default_project_name = "ribosome profiling data"
default_min_metagene_profile_count = 1000
default_min_metagene_profile_bayes_factor_mean = 5
default_max_metagene_profile_bayes_factor_var = 5
default_min_visualization_count = 500
default_num_cpus = 1

abstract = """
This document shows the results of preprocessing steps. 
In particular, it shows the amount of reads filtered at each step in the pipeline. 
Additionally, it includes ``metagene'' profile figures for all of the samples.
"""

commands = {
    "riboseq": r"\textsc{ribo-seq}\xspace"
}

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
Figure~\\ref{fig:mapping-info}(right) shows a ``zoomed in'' version which does not include the reads of poor quality and that mapped to ribosomal sequences.
"""

read_filtering_caption = """
The number of reads passing each phase of the filtering pipeline.
"""

def create_fastqc_reports(name_data, config, args):
    name, data = name_data
    msg = "{}: creating fastqc reports".format(name)
    logger.info(msg)

    note = config.get('note', None)

    # keep multimappers?
    is_unique = not ('keep_riboseq_multimappers' in config)

    # first, get the ribo_filenames
    raw_data = data
    without_adapters = filenames.get_without_adapters_fastq(
        config['riboseq_data'], name, note=note)
    with_rrna = filenames.get_with_rrna_fastq(
        config['riboseq_data'], name, note=note)
    without_rrna = filenames.get_without_rrna_fastq(
        config['riboseq_data'], name, note=note)
    genome_bam = filenames.get_riboseq_bam(
        config['riboseq_data'], name, note=note)
    unique_bam = filenames.get_riboseq_bam(
        config['riboseq_data'], name, is_unique=is_unique, note=note)

    if args.create_fastqc_reports:

        # now, get the fastqc report ribo_filenames
        raw_data_fastqc = filenames.get_raw_data_fastqc_data(
            config['riboseq_data'], raw_data)
        without_adapters_fastqc = filenames.get_without_adapters_fastqc_data(
            config['riboseq_data'], name, note=note)
        with_rrna_fastqc = filenames.get_with_rrna_fastqc_data(
            config['riboseq_data'], name, note=note)
        without_rrna_fastqc = filenames.get_without_rrna_fastqc_data(
            config['riboseq_data'], name, note=note)

        genome_bam_fastqc = filenames.get_riboseq_bam_fastqc_data(
            config['riboseq_data'], name, note=note)
        unique_bam_fastqc = filenames.get_riboseq_bam_fastqc_data(
            config['riboseq_data'], name, is_unique=is_unique, note=note)

        # create the fastqc reports if they do not exist
        raw_data_fastqc_path = filenames.get_raw_data_fastqc_path(config['riboseq_data'])
        without_adapters_fastqc_path = filenames.get_without_adapters_fastqc(config['riboseq_data'])
        with_rrna_fastqc_path = filenames.get_with_rrna_fastqc(config['riboseq_data'])
        without_rrna_fastqc_path = filenames.get_without_rrna_fastqc(config['riboseq_data'])
        without_rrna_mapping_fastqc_path = filenames.get_riboseq_bam_fastqc_path(config['riboseq_data'])

        fastqc_tmp_str = ""
        if args.tmp is not None:
            fastqc_tmp_str = "--dir {}".format(args.tmp)

        msg = "Looking for raw data fastqc report: '{}'".format(raw_data_fastqc)
        logger.debug(msg)
        cmd = "fastqc --outdir {} --extract {} {}".format(raw_data_fastqc_path, 
            raw_data, fastqc_tmp_str)
        in_files = [raw_data]
        out_files = [raw_data_fastqc]

        shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite)

        cmd = "fastqc --outdir {} --extract {} {}".format(without_adapters_fastqc_path, 
            without_adapters, fastqc_tmp_str)
        in_files = [without_adapters]
        out_files = [without_adapters_fastqc]
        shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite)

        cmd = "fastqc --outdir {} --extract {} {}".format(with_rrna_fastqc_path, 
            with_rrna, fastqc_tmp_str)
        in_files = [with_rrna]
        out_files = [with_rrna_fastqc]
        shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite)

        cmd = "fastqc --outdir {} --extract {} {}".format(without_rrna_fastqc_path, 
            without_rrna, fastqc_tmp_str)
        in_files = [without_rrna]
        out_files = [without_rrna_fastqc]
        shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite)

        cmd = "fastqc --outdir {} --extract {} {}".format(without_rrna_mapping_fastqc_path, 
            genome_bam, fastqc_tmp_str)
        in_files = [genome_bam]
        out_files = [genome_bam_fastqc]

        msg = "genome_bam: '{}'".format(genome_bam)
        logger.debug(msg)
        
        msg = "genome_bam_fastqc: '{}'".format(genome_bam_fastqc)
        logger.debug(msg)

        shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite)

        cmd = "fastqc --outdir {} --extract {} {}".format(without_rrna_mapping_fastqc_path, 
            unique_bam, fastqc_tmp_str)
        in_files = [unique_bam]
        out_files = [unique_bam_fastqc]
        shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files, overwrite=args.overwrite)

        # in some cases, fastqc can fail. make sure all of the reports are present
        all_fastqc_reports = [
            raw_data_fastqc,
            without_adapters_fastqc,
            without_rrna_fastqc,
            genome_bam_fastqc,
            unique_bam_fastqc
        ]

        missing_files = [
            f for f in all_fastqc_reports if not os.path.exists(f)
        ]

        if len(missing_files) > 0:
            msg = "The following fastqc reports were not created correctly:\n"
            msg += '\n'.join(missing_files)
            logger.warning(msg)


def create_figures(config_file, config, name, offsets_df, args):
    """ This function creates all of the figures in the preprocessing report
        for the given dataset.
    """
    logging_str = logging_utils.get_logging_options_string(args)
    note = config.get('note', None)
    
    note_str = filenames.get_note_string(note)

    # keep multimappers?
    is_unique = not ('keep_riboseq_multimappers' in config)

    image_type_str = "--image-type {}".format(args.image_type)

    min_read_length = int(offsets_df['length'].min())
    max_read_length = int(offsets_df['length'].max())

    min_read_length_str = "--min-read-length {}".format(min_read_length)
    max_read_length_str = "--max-read-length {}".format(max_read_length)

    msg = "{}: Getting and visualizing read length distribution".format(name)
    logger.info(msg)

    # all aligned reads
    genome_bam = filenames.get_riboseq_bam(
        config['riboseq_data'], name, note=note)
            
    # uniquely aligned reads
    unique_filename = filenames.get_riboseq_bam(
        config['riboseq_data'], name, is_unique=is_unique, note=note)

    # the read length counts
    read_length_distribution = filenames.get_riboseq_read_length_distribution(
        config['riboseq_data'], name, note=note)
    
    # the plots
    cmd = "get-read-length-distribution {} {} --out {} {}".format(
        genome_bam, 
        unique_filename,
        read_length_distribution, 
        logging_str
    )
    in_files = [genome_bam, unique_filename]
    out_files = [read_length_distribution]
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files, 
        overwrite=args.overwrite, call=True)

    # visualize all read counts
    title = None
    if 'riboseq_sample_name_map' in config:
        title = config['riboseq_sample_name_map'].get(name)
    if title is None:
        title = "{}{}".format(name, note_str)
    
    title_str = "{}, All aligned reads".format(title)        
    title_str = "--title={}".format(shlex.quote(title_str))

    # get the basename for the distribution file
    unique_str = filenames.get_unique_string(False)
    sample_name = "{}{}{}".format(name, note_str, unique_str)

    read_length_distribution_image = filenames.get_riboseq_read_length_distribution_image(
        config['riboseq_data'], 
        name, 
        is_unique=False, 
        note=note, 
        image_type=args.image_type
    )
    
    cmd = "plot-read-length-distribution {} {} {} {} {} {}".format(
        read_length_distribution, 
        sample_name,
        read_length_distribution_image, 
        title_str, 
        min_read_length_str, 
        max_read_length_str
    )

    in_files = [read_length_distribution]
    out_files = [read_length_distribution_image]
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files, 
        overwrite=args.overwrite, call=True)
    
    # visualize unique read counts

    # we already have the title
    title_str = "{}, Uniquely aligned reads".format(title)        
    title_str = "--title={}".format(shlex.quote(title_str))

    unique_read_length_distribution_image = filenames.get_riboseq_read_length_distribution_image(
        config['riboseq_data'], 
        name, 
        is_unique=is_unique, 
        note=note, 
        image_type=args.image_type
    )

    # get the basename for the distribution file
    unique_str = filenames.get_unique_string(True)
    sample_name = "{}{}{}".format(name, note_str, unique_str)

    cmd = "plot-read-length-distribution {} {} {} {} {} {}".format(
        read_length_distribution, 
        sample_name,
        unique_read_length_distribution_image, 
        title_str, 
        min_read_length_str, 
        max_read_length_str
    )
    in_files = [read_length_distribution]
    out_files = [unique_read_length_distribution_image]
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files, 
        overwrite=args.overwrite, call=True)

    # visualize the metagene profiles
    msg = "{}: Visualizing metagene profiles and Bayes' factors".format(name)
    logger.info(msg)

    metagene_profiles = filenames.get_metagene_profiles(config['riboseq_data'], 
        name, is_unique=is_unique, note=note)
    
    profile_bayes_factor = filenames.get_metagene_profiles_bayes_factors(
        config['riboseq_data'],
        name, 
        is_unique=is_unique, 
        note=note
    )

    mp_df = pd.read_csv(metagene_profiles)

    for length in range(min_read_length, max_read_length+1):

        mask_length = offsets_df['length'] == length

        # make sure we had some reads of that length
        if sum(mask_length) == 0:
            continue
        length_row = offsets_df[mask_length].iloc[0]
               
        # make sure we have enough reads to visualize
        if length_row['highest_peak_profile_sum'] < args.min_visualization_count:
            continue
        
        # visualize the metagene profile
        metagene_profile_image = filenames.get_metagene_profile_image(
            config['riboseq_data'], 
            name, 
            image_type=args.image_type, 
            is_unique=is_unique, 
            length=length, 
            note=note
        )

        title_str = "{}. length: {}".format(title, length)
        title_str = "--title {}".format(shlex.quote(title_str))
        cmd = ("create-read-length-metagene-profile-plot {} {} {} {}".format(
            metagene_profiles, length, metagene_profile_image, title_str))
        in_files = [metagene_profiles]
        out_files = [metagene_profile_image]
        shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files, 
            overwrite=args.overwrite, call=True)

        # and the Bayes' factor
        if args.show_read_length_bfs:
            metagene_profile_image = filenames.get_metagene_profile_bayes_factor_image(
                config['riboseq_data'], 
                name, 
                image_type=args.image_type, 
                is_unique=is_unique, 
                length=length, 
                note=note
            )

            title_str = "Metagene profile Bayes' factors: {}. length: {}".format(title, length)
            title_str = "--title {}".format(shlex.quote(title_str))
            fontsize_str = "--font-size 15"

            cmd = ("visualize-metagene-profile-bayes-factor {} {} {} {} {}".format(
                profile_bayes_factor, 
                length, 
                metagene_profile_image, 
                title_str,
                fontsize_str
            ))

            in_files = [profile_bayes_factor]
            out_files = [metagene_profile_image]
            shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files, 
                overwrite=args.overwrite, call=True)

    # the orf-type metagene profiles
    if args.show_orf_periodicity:
        msg = "{}: Visualizing the ORF type metagene profiles".format(title)
        logger.info(msg)


        try:
            lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(config, 
                name, is_unique=is_unique)
        except FileNotFoundError:
            msg = ("Could not parse out lengths and offsets for sample: {}. "
                "Skipping".format(name))
            logger.error(msg)
            return


        orfs_genomic = filenames.get_orfs(config['genome_base_path'], 
            config['genome_name'], note=config.get('orf_note'))
             
        profiles = filenames.get_riboseq_profiles(config['riboseq_data'], name, 
                length=lengths, offset=offsets, is_unique=is_unique, note=note)

        title_str = "{}, ORF-type periodicity".format(title)
        title_str = "--title {}".format(shlex.quote(title_str))
        
        orf_type_profile_base = filenames.get_orf_type_profile_base(
            config['riboseq_data'], name, length=lengths, offset=offsets, 
            is_unique=is_unique, note=note, subfolder='orf-profiles')

        strand = "+"
        orf_type_profiles_forward = [
            filenames.get_orf_type_profile_image(
                orf_type_profile_base, 
                orf_type, 
                strand, 
                args.image_type
            )   for orf_type in ribo_utils.orf_types
        ]
        
        strand = "-"
        orf_type_profiles_reverse = [
            filenames.get_orf_type_profile_image(
                orf_type_profile_base, 
                orf_type, 
                strand, 
                args.image_type
            )   for orf_type in ribo_utils.orf_types
        ]

        cmd = ("visualize-orf-type-metagene-profiles {} {} {} {} {} {}".format(
            orfs_genomic, profiles, orf_type_profile_base, title_str, 
            image_type_str, logging_str))

        in_files = [orfs_genomic, profiles]
        out_files = orf_type_profiles_forward + orf_type_profiles_reverse
        shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files, 
            overwrite=args.overwrite)


def create_read_filtering_plots(config_file, config, args):
        
    # get the filtering counts
    note = config.get('note', None)
    read_filtering_counts = filenames.get_riboseq_read_filtering_counts(
        config['riboseq_data'], note=note)
    overwrite_str = ""
    if args.overwrite:
        overwrite_str = "--overwrite"

    logging_str = logging_utils.get_logging_options_string(args)

    cpus_str = "--num-cpus {}".format(args.num_cpus)
    cmd = "get-all-read-filtering-counts {} {} {} {} {}".format(
        config_file, 
        read_filtering_counts, 
        overwrite_str, 
        cpus_str, 
        logging_str
    )

    in_files = [config_file]
    out_files = [read_filtering_counts]
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files, 
        overwrite=args.overwrite, call=True)

    # and visualize them
    read_filtering_image = filenames.get_riboseq_read_filtering_counts_image(
        config['riboseq_data'], note=note, image_type=args.image_type)
    
    title = "Read filtering counts"
    title_str = "--title {}".format(shlex.quote(title))
    cmd = "visualize-read-filtering-counts {} {} {} --config {}".format(
        read_filtering_counts, 
        read_filtering_image, 
        title_str,
        config_file
    )
    in_files = [read_filtering_counts]
    out_files=[read_filtering_image]
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files, 
        overwrite=args.overwrite, call=True)

    # and visualize the filtering without the rrna
    n = "no-rrna-{}".format(note)
    read_filtering_image = filenames.get_riboseq_read_filtering_counts_image(
        config['riboseq_data'], note=n, image_type=args.image_type)
    
    title = "Read filtering counts, no ribosomal matches"
    title_str = "--title {}".format(shlex.quote(title))
    cmd = "visualize-read-filtering-counts {} {} {} --config {} --without-rrna".format(
        read_filtering_counts, 
        read_filtering_image, 
        title_str,
        config_file
    )

    in_files = [read_filtering_counts]
    out_files=[read_filtering_image]
    shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files, 
        overwrite=args.overwrite, call=True)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script creates a simple latex document containing the read "
        "filtering images, metagene profiles and analysis, and standard section text.")
    parser.add_argument('config', help="The (yaml) config file for the project")
    parser.add_argument('out', help="The path for the output files")

    parser.add_argument('--show-orf-periodicity', help="If this flag is "
        "present, bar charts showing the periodicity of each ORF type will be "
        "included in the report.", action='store_true')

    parser.add_argument('--show-read-length-bfs', help="If this flag is given, "
        "plots showing the Bayes factor at each offset for each read length "
        "are included in the report.", action='store_true')

    parser.add_argument('--overwrite', help="If this flag is present, existing files will "
        "be overwritten.", action='store_true')

    parser.add_argument('--min-visualization-count', help="Read lengths with fewer than this "
        "number of reads will not be included in the report.", type=int, 
        default=default_min_visualization_count)

    parser.add_argument('--image-type', help="The type of image types to create. This "
        "must be an extension which matplotlib can interpret.", default=default_image_type)

    parser.add_argument('-c', '--create-fastqc-reports', help="If this flag is given, then "
        "fastqc reports will be created for most fastq and bam files. By default, they are "
        "not created.", action='store_true')
    
    parser.add_argument('--tmp', help="If the fastqc reports are created, "
        "they will use this location for temp files", default=default_tmp)
     
    parser.add_argument('--note', help="If this option is given, it will be used in the "
        "filenames.\n\nN.B. This REPLACES the note in the config file.", default=default_note)

    slurm.add_sbatch_options(parser)
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    config = yaml.load(open(args.config))

    if args.note is not default_note:
        config['note'] = args.note
    note = config.get('note', None)
    
    sample_names = sorted(config['riboseq_samples'].keys())
    
    # keep multimappers?
    is_unique = not ('keep_riboseq_multimappers' in config)

    programs =  [   
	'create-read-length-metagene-profile-plot',
        'visualize-metagene-profile-bayes-factor',
        'get-all-read-filtering-counts',
        'samtools',
        'visualize-read-filtering-counts',
        'get-read-length-distribution',
        'plot-read-length-distribution'
    ]

    if args.create_fastqc_reports:
        programs.extend(['fastqc','java'])
        
    shell_utils.check_programs_exist(programs)

    if args.use_slurm:
        cmd = ' '.join(sys.argv)
        slurm.check_sbatch(cmd, args=args)
        return

    # make sure the path to the output file exists
    os.makedirs(args.out, exist_ok=True)
   
    # first, create the read filtering information... 
    create_read_filtering_plots(args.config, config, args)
    # ... and all the other figures.
    for name in sample_names:
        periodic_offsets = filenames.get_periodic_offsets(config['riboseq_data'], 
            name, is_unique=is_unique, note=note)
        offsets_df = pd.read_csv(periodic_offsets)
        create_figures(args.config, config, name, offsets_df, args)

    min_metagene_profile_count = config.get(
        "min_metagene_profile_count", default_min_metagene_profile_count)

    min_metagene_profile_bayes_factor_mean = config.get(
        "min_metagene_profile_bayes_factor_mean", 
        default_min_metagene_profile_bayes_factor_mean)

    max_metagene_profile_bayes_factor_var = config.get(
        "max_metagene_profile_bayes_factor_var", 
        default_max_metagene_profile_bayes_factor_var)

    project_name = config.get("project_name", default_project_name)
    title = "Preprocessing results for {}".format(project_name)

    tex_file = os.path.join(args.out, "preprocessing-report.tex")
    with open(tex_file, 'w') as out:

        latex.begin_document(out, title, abstract, commands=commands)

        latex.section(out, "Introduction")

        latex.clearpage(out)
        latex.newpage(out)

        latex.section(out, "Mapping and filtering")
        latex.write(out, mapping_and_filtering_text)

        # the read filtering figures
        read_filtering_image = filenames.get_riboseq_read_filtering_counts_image(
            config['riboseq_data'], note=note, image_type=args.image_type)
    
        n = "no-rrna-{}".format(note)
        no_rrna_read_filtering_image = filenames.get_riboseq_read_filtering_counts_image(
            config['riboseq_data'], note=n, image_type=args.image_type)

        latex.begin_figure(out)
        latex.write_graphics(out, read_filtering_image, width=0.45)
        latex.write_graphics(out, no_rrna_read_filtering_image, width=0.45)
        latex.write_caption(out, read_filtering_caption, label=read_filtering_label)
        latex.end_figure(out)

        latex.clearpage(out)

        # the read length distributions
        latex.section(out, "Read length distributions", 
            label=length_distribution_section_label)

        msg = "Writing length distribution figures"
        logger.info(msg)

        latex.begin_table(out, "cc")

        latex.write_header(out, ["All aligned reads", "Uniquely-aligning reads"])

        for name in sample_names:
            data = config['riboseq_samples'][name]
            read_length_distribution_image = filenames.get_riboseq_read_length_distribution_image(
                config['riboseq_data'], name, is_unique=False, note=note, 
                image_type=args.image_type)

            unique_read_length_distribution_image = filenames.get_riboseq_read_length_distribution_image(
                config['riboseq_data'], name, is_unique=True, note=note, 
                image_type=args.image_type)

            
            msg = "Looking for image file: {}".format(read_length_distribution_image)
            logger.debug(msg)

            if os.path.exists(read_length_distribution_image):
                latex.write_graphics(out, read_length_distribution_image, width=0.45)
            else:
                msg = "Could not find image: {}".format(read_length_distribution_image)
                logger.warning(msg)
                
                text = "Missing: {}\n\n".format(name)
                latex.write(out, text)

            latex.write_column_sep(out)
            

            msg = "Looking for image file: {}".format(unique_read_length_distribution_image)
            logger.debug(msg)

            if os.path.exists(unique_read_length_distribution_image):
                latex.write_graphics(out, unique_read_length_distribution_image, width=0.45)
            else:
                msg = "Could not find image: {}".format(unique_read_length_distribution_image)
                logger.warning(msg)
            
                text = "Missing: {}\n\n".format(name)
                latex.write(out, text)

            latex.write_row_sep(out)

            

        latex.end_table(out)
        latex.clearpage(out)


        latex.section(out, "Read length periodicity", label=periodicity_label)

        for name in sample_names:
            i = 0

            data = config['riboseq_samples'][name]

            msg = "Processing sample: {}".format(name)
            logger.info(msg)

            logger.debug("overwrite: {}".format(args.overwrite))
    
            periodic_offsets = filenames.get_periodic_offsets(config['riboseq_data'], 
                name, is_unique=is_unique, note=note)
            offsets_df = pd.read_csv(periodic_offsets)

            min_read_length = int(offsets_df['length'].min())
            max_read_length = int(offsets_df['length'].max())
    

            latex.begin_table(out, "YY")

            header = "\\multicolumn{2}{c}{" + name + "}"
            header = [header]
            latex.write_header(out, header)

            for length in range(min_read_length, max_read_length + 1):
                msg = "Processing length: {}".format(length)
                logger.info(msg)

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
                    logger.warning(msg)
                    continue

                metagene_profile_image = filenames.get_metagene_profile_image(
                    config['riboseq_data'], name, image_type=args.image_type, 
                    is_unique=is_unique, length=length, note=note)
                
                #title = ("length: {}. P-site offset: {}. \\newline status: {}"
                    #"\n".format(length, offset, offset_status))
                #latex.write(out, title, size="scriptsize")
                title = ("Length: {}. P-site offset: {}. Status: {}\n".format(length, offset, offset_status))
                if args.show_read_length_bfs:
                    title = "\scriptsize{" + title + "}"
                    title = "\\multicolumn{2}{c}{" + title + "}"
                    latex.write(out, title)
                    latex.write_row_sep(out)
                else:
                    latex.write(out, title, size="scriptsize")

                latex.write_graphics(out, metagene_profile_image, width=0.45)
                               
                i += 1
                if i%2 == 1:
                    latex.write_column_sep(out)
                else:
                    latex.write_row_sep(out)


                if args.show_read_length_bfs:
                    
                    bayes_factor_image = filenames.get_metagene_profile_bayes_factor_image(
                        config['riboseq_data'], name, image_type=args.image_type, 
                        is_unique=is_unique, length=length, note=note)

                    #latex.centering(out)
                    latex.write_graphics(out, bayes_factor_image, width=0.45)
                        
                    i += 1
                    if i%2 == 1:
                        latex.write_column_sep(out)
                    else:
                        latex.write_row_sep(out)

                               
            if i%2 == 1:
                latex.write_row_sep(out)

            latex.end_table(out)
            latex.clearpage(out)

        ### ORF type metagene profiles
        if args.show_orf_periodicity:
            title = "ORF type periodicity"
            latex.section(out, title)
            
            strands = ['+', '-']
            for sample_name in sample_names:
                i = 0

                try:
                    lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(
                        config, sample_name, is_unique=is_unique)
                except FileNotFoundError:
                    msg = ("Could not parse out lengths and offsets for sample: {}. "
                        "Skipping".format(sample_name))
                    logger.error(msg)
                    continue
                
                orf_type_profile_base = filenames.get_orf_type_profile_base(
                    config['riboseq_data'], sample_name, length=lengths, offset=offsets, 
                    is_unique=is_unique, note=note, subfolder='orf-profiles')

                for orf_type in ribo_utils.orf_types:
                    for strand in strands:
                        orf_type_profile = filenames.get_orf_type_profile_image(
                            orf_type_profile_base, orf_type, strand, 
                            image_type=args.image_type)


                        msg = "Looking for image file: {}".format(orf_type_profile)
                        logger.debug(msg)
                        if os.path.exists(orf_type_profile):
                            if i%4 == 0:
                                latex.begin_figure(out)

                            i += 1
                            latex.write_graphics(out, orf_type_profile, height=0.23)

                            if i%4 == 0:
                                latex.end_figure(out)
                                latex.clearpage(out)

                if (i>0) and (i%4 != 0):
                    latex.end_figure(out)
                    latex.clearpage(out)

        latex.end_document(out)

    tex_filename = os.path.basename(tex_file)
    latex.compile(args.out, tex_filename)

    if args.create_fastqc_reports:
        parallel.apply_parallel_iter(config['riboseq_samples'].items(), 
            args.num_cpus, 
            create_fastqc_reports, config, args)


if __name__ == '__main__':
    main()

