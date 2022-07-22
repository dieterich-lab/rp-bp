#! /usr/bin/env python3

import argparse
import itertools
import logging
import os
import shlex
import sys
import yaml

import pbio.misc.latex as latex
import pbio.misc.logging_utils as logging_utils
import pbio.misc.parallel as parallel
import pbio.misc.shell_utils as shell_utils
import pbio.misc.slurm as slurm
import pbio.misc.utils as utils

import pbio.ribo.ribo_filenames as filenames
import pbio.ribo.ribo_utils as ribo_utils

from rpbp.defaults import default_num_cpus, metagene_options

logger = logging.getLogger(__name__)

default_project_name = "riboseq"
default_image_type = 'eps'

default_uniprot = ""
default_uniprot_label = "UniRef"


def _create_figures(name_pretty_name_is_replicate, config, args):
    """ This function creates all of the figures in the prediction report
        for the given dataset.
    """
    name, pretty_name, is_replicate = name_pretty_name_is_replicate
    
    # keep multimappers?
    is_unique = not ('keep_riboseq_multimappers' in config)

    # by default, we will not include chisq
    chisq_values = [False]
    if args.show_chisq:
        chisq_values = [True, False]

    filtered_values = [True]
    if args.show_unfiltered_orfs:
        filtered_values = [True, False]

    grouped_values = [True, False]

    logging_str = logging_utils.get_logging_options_string(args)

    note_str = config.get('note', None)
    out_note_str = config.get('note', None)
    if args.note is not None and len(args.note) > 0:
        out_note_str = args.note

    image_type_str = "--image-type {}".format(args.image_type)
    num_cpus_str = "--num-cpus {}".format(args.num_cpus)
   
    fraction = config.get('smoothing_fraction', None)
    reweighting_iterations = config.get('smoothing_reweighting_iterations', None)

    # if this is a replicate, we do not worry about lengths and offsets
    if is_replicate:
        lengths = None
        offsets = None
    else:
        try:
            lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(config,
                                                                           name,
                                                                           is_unique=is_unique,
                                                                           default_params=metagene_options)
        except FileNotFoundError:
            msg = ("Could not parse out lengths and offsets for sample: {}. "
                "Skipping".format(name))
            logger.error(msg)
            return
        
    unsmoothed_profiles = filenames.get_riboseq_profiles(
        config['riboseq_data'],
        name,
        length=lengths,
        offset=offsets,
        is_unique=is_unique,
        note=note_str,
        is_smooth=False
    )

    msg = "{}: creating the ORF types bar charts".format(name)
    logger.debug(msg)

    it = itertools.product(grouped_values, chisq_values, filtered_values)

    for is_grouped, is_chisq, is_filtered in it:

        is_grouped_str = ""
        if is_grouped:
            is_grouped_str = ", Grouped"

        is_filtered_str = ""
        if is_filtered:
            is_filtered_str = ", Filtered"
        
        if is_chisq:
            title_str = "{}{}{}, Rp-$\chi^2$".format(pretty_name, is_grouped_str, is_filtered_str)
            title_str = shlex.quote(title_str)
            title_str = "--title {}".format(title_str)

            f = None
            rw = None

            orfs = filenames.get_riboseq_predicted_orfs(config['riboseq_data'], 
                name, length=lengths, offset=offsets, is_unique=is_unique, note=note_str, 
                is_chisq=True, is_filtered=is_filtered)

        else:
            title_str = "{}{}{}, Rp-Bp".format(pretty_name, is_grouped_str, is_filtered_str)
            title_str = shlex.quote(title_str)
            title_str = "--title {}".format(title_str)

            f = fraction
            rw = reweighting_iterations
            orfs = filenames.get_riboseq_predicted_orfs(config['riboseq_data'], 
                name, length=lengths, offset=offsets, is_unique=is_unique, note=note_str, 
                fraction=f, reweighting_iterations=rw, 
                is_filtered=is_filtered)

        use_groups_str = ""
        if is_grouped:
            use_groups_str = "--use-groups"
        
        orf_types_bar_chart = filenames.get_orf_types_bar_chart(
            config['riboseq_data'], 
            name, 
            length=lengths, 
            offset=offsets, 
            is_unique=is_unique, 
            note=out_note_str, 
            image_type=args.image_type,
            fraction=f, 
            reweighting_iterations=rw,
            is_grouped=is_grouped, 
            is_chisq=is_chisq, 
            is_filtered=is_filtered
        )

        cmd = "create-orf-types-bar-chart {} {} {} {}".format(
            orfs, 
            orf_types_bar_chart,
            title_str, 
            use_groups_str
        )

        in_files = [orfs]
        out_files = [orf_types_bar_chart]
        shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files, 
            overwrite=args.overwrite)

    
    msg = "{}: creating the ORF length distributions line graph".format(name)
    logger.debug(msg)

    uniprot_str = ""
    uniprot_label_str = ""
    if os.path.exists(args.uniprot):
        uniprot_str = "--uniprot {}".format(args.uniprot)
        uniprot_label_str = shlex.quote(args.uniprot_label)
        uniprot_label_str = "--uniprot-label {}".format(uniprot_label_str)

    for is_grouped in grouped_values:
        for is_chisq in chisq_values:
            
            if is_chisq:
                title_str = "{}, Rp-$\chi^2$".format(pretty_name)
                title_str = shlex.quote(title_str)
                title_str = "--title {}".format(title_str)

                f = None
                rw = None
                
                orfs = filenames.get_riboseq_predicted_orfs(
                    config['riboseq_data'], 
                    name, 
                    length=lengths, 
                    offset=offsets, 
                    is_unique=is_unique, 
                    note=note_str, 
                    is_chisq=True
                )

            else:
                title_str = "{}, Rp-Bp".format(pretty_name)
                title_str = shlex.quote(title_str)
                title_str = "--title {}".format(title_str)

                f = fraction
                rw = reweighting_iterations

                orfs = filenames.get_riboseq_predicted_orfs(
                    config['riboseq_data'], 
                    name, 
                    length=lengths, 
                    offset=offsets, 
                    is_unique=is_unique, 
                    note=note_str, 
                    fraction=f, 
                    reweighting_iterations=rw
                )


            use_groups_str = ""
            if is_grouped:
                use_groups_str = "--use-groups"
            
            orf_length_line_graph = filenames.get_orf_length_distribution_line_graph(
                config['riboseq_data'], 
                name, 
                length=lengths, 
                offset=offsets, 
                is_unique=is_unique, 
                note=out_note_str, 
                image_type=args.image_type,
                fraction=f, 
                reweighting_iterations=rw,
                is_grouped=is_grouped, 
                is_chisq=is_chisq
            )

            cmd = ("create-orf-length-distribution-line-graph {} {} {} {} {} {}".format(
                orfs, 
                orf_length_line_graph, 
                title_str, 
                use_groups_str, 
                uniprot_str, 
                uniprot_label_str
            ))

            in_files = [orfs]
            out_files = [orf_length_line_graph]
            shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files, 
                overwrite=args.overwrite)

    if args.show_orf_periodicity:
        msg = "{}: creating the ORF type metagene profiles".format(name)
        logger.debug(msg)

        for is_chisq in chisq_values:
            
            if is_chisq:
                title_str = "{}, Rp-$\chi^2$".format(pretty_name)
                title_str = shlex.quote(title_str)
                title_str = "--title {}".format(title_str)
                f = None
                rw = None
                is_smooth = False
                profiles = unsmoothed_profiles
        
                orfs = filenames.get_riboseq_predicted_orfs(
                    config['riboseq_data'], 
                    name, 
                    length=lengths, 
                    offset=offsets, 
                    is_unique=is_unique, 
                    note=note_str, 
                    is_chisq=True, 
                    is_filtered=is_filtered
                )

            else:
                title_str = "{}, Rp-Bp".format(pretty_name)
                title_str = shlex.quote(title_str)
                title_str = "--title {}".format(title_str)

                f = fraction
                rw = reweighting_iterations
                is_smooth = False
                profiles = unsmoothed_profiles

                orfs = filenames.get_riboseq_predicted_orfs(
                    config['riboseq_data'], 
                    name, 
                    length=lengths, 
                    offset=offsets, 
                    is_unique=is_unique, 
                    note=note_str, 
                    fraction=f, 
                    reweighting_iterations=rw
                )

            
            orf_type_profile_base = filenames.get_orf_type_profile_base(
                config['riboseq_data'], 
                name, 
                length=lengths, 
                offset=offsets, 
                is_unique=is_unique, 
                note=out_note_str, 
                fraction=f, 
                reweighting_iterations=rw,
                is_chisq=is_chisq
            )

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
                orfs,
                profiles,
                orf_type_profile_base,
                title_str, 
                image_type_str,
                logging_str
            ))

            in_files = [orfs]
            out_files = orf_type_profiles_forward + orf_type_profiles_reverse
            shell_utils.call_if_not_exists(cmd, out_files, in_files=in_files, 
                overwrite=args.overwrite)

def create_all_figures(config, args):
    
    # first, handle all of the regular datasets

    is_replicate = False
    sample_names = sorted(config['riboseq_samples'].keys())
    sample_name_map = ribo_utils.get_sample_name_map(config)
    samples = [(name, sample_name_map[name], is_replicate) 
        for name in sample_names
    ]
    
    is_replicate = True
    replicate_names = sorted(ribo_utils.get_riboseq_replicates(config).keys())
    condition_name_map = ribo_utils.get_riboseq_condition_name_map(config)
    conditions = [(name, condition_name_map[name], is_replicate)
        for name in replicate_names
    ]

    all_names = samples + conditions
    parallel.apply_parallel_iter(
        all_names,
        args.num_cpus,
        _create_figures,
        config,
        args,
        progress_bar=True
    )

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script creates the plots which detail the basic characteristics "
        "of the ORF predictions from the Rp-Bp pipeline. It also creates and compiles (if "
        "possible) a latex report for them.")

    parser.add_argument('config', help="The (yaml) config file")
    parser.add_argument('out', help="The base output directory for the latex report")

    parser.add_argument('--show-unfiltered-orfs', help="If this flag is "
        "present, bar charts showing the distribution of the types of the "
        "unfiltered ORF set will be included", action='store_true')
    
    parser.add_argument('--show-orf-periodicity', help="If this flag is "
        "present, bar charts showing the periodicity of each ORF type will be "
        "included in the report.", action='store_true')

    parser.add_argument('--uniprot', help="The uniprot ORF lengths, if available", 
        default=default_uniprot)
    parser.add_argument('--uniprot-label', help="The label to use for the uniprot ORFs in "
        "the plot", default=default_uniprot_label)
    
    parser.add_argument('--image-type', help="The format of the image files. This must be "
        "a format usable by matplotlib.", default=default_image_type)

    parser.add_argument('--overwrite', help="If this flag is present, existing files will "
        "be overwritten.", action='store_true')
        
    parser.add_argument('--note', help="If this option is given, it will be used in the "
        "filenames.\n\nN.B. This REPLACES the note in the config file.", default=None)

    parser.add_argument('--show-chisq', help="If this flag is given, then the "
        "results from Rp-chi will be included in the document; otherwise, they "
        "will not be created or shown.", action='store_true')

    parser.add_argument('-t', '--tmp', help="A location for temporary files",
        default=None)

    slurm.add_sbatch_options(parser, num_cpus=default_num_cpus)
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    config = yaml.load(open(args.config), Loader=yaml.FullLoader)
    
    # keep multimappers?
    is_unique = not ('keep_riboseq_multimappers' in config)

    programs =  [   
        'create-orf-length-distribution-line-graph',
        'create-orf-types-bar-chart',
        'visualize-orf-type-metagene-profiles'
    ]
    shell_utils.check_programs_exist(programs)
    
    required_keys = [
        'riboseq_data',
        'riboseq_samples'
    ]
    utils.check_keys_exist(config, required_keys)

    if args.use_slurm:
        cmd = ' '.join(sys.argv)
        slurm.check_sbatch(cmd, args=args)
        return

    # by default, we will not include chisq
    chisq_values = [False]
    if args.show_chisq:
        chisq_values = [True, False]

    filtered_values = [True]
    if args.show_unfiltered_orfs:
        filtered_values = [True, False]

    grouped_values = [True, False]
    
    # make sure the path to the output file exists
    os.makedirs(args.out, exist_ok=True)

    # first, create all of the figures
    create_all_figures(config, args)

    note_str = config.get('note', None)
    out_note_str = note_str

    if args.note is not None and len(args.note) > 0:
        out_note_str = args.note
    
    fraction = config.get('smoothing_fraction', None)
    reweighting_iterations = config.get('smoothing_reweighting_iterations', None)

    project_name = config.get("project_name", default_project_name)
    title = "Rp-Bp prediction analysis for {}".format(project_name)
    abstract = "This document shows the results of the Rp-Bp pipeline analysis."

    
    #tex_file = os.path.join(args.out, "prediction-report.tex")
    tex_file = filenames.get_rpbp_prediction_report(args.out, out_note_str)

    
    with open(tex_file, 'w') as out:

        latex.begin_document(out, title, abstract)

        latex.write(out, "\n")

        latex.clearpage(out)

        ### ORF type distributions
        title = "Predicted ORF type distributions"
        latex.section(out, title)

        # first, handle all of the regular datasets
        sample_names = sorted(config['riboseq_samples'].keys())
        
        # and check if we also have replicates
        replicate_names = []
        if 'riboseq_biological_replicates' in config:
            replicate_names = sorted(ribo_utils.get_riboseq_replicates(config).keys())

        strands = ["+", "-"]

        i = 0
        for sample_name in sample_names:
            
            try:
                lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(config,
                                                                               sample_name,
                                                                               is_unique=is_unique,
                                                                               default_params=metagene_options)
            except FileNotFoundError:
                msg = ("Could not parse out lengths and offsets for sample: {}. "
                    "Skipping".format(sample_name))
                logger.error(msg)
                continue
            
            caption = "ORF types: {}".format(sample_name)
            is_first = True

            # first, just dump all of the bar charts to the page            
            it = itertools.product(grouped_values, chisq_values, filtered_values)

            for is_grouped, is_chisq, is_filtered in it:

                if is_chisq:
                    f = None
                    rw = None
                else:
                    f = fraction
                    rw = reweighting_iterations
                
                orf_types_bar_chart = filenames.get_orf_types_bar_chart(
                    config['riboseq_data'], 
                    sample_name, 
                    length=lengths, 
                    offset=offsets, 
                    is_unique=is_unique, 
                    note=out_note_str, 
                    image_type=args.image_type,
                    fraction=f, 
                    reweighting_iterations=rw,
                    is_grouped=is_grouped, 
                    is_chisq=is_chisq, 
                    is_filtered=is_filtered
                )

                msg = "Looking for image file: {}".format(orf_types_bar_chart)
                logger.debug(msg)

                if os.path.exists(orf_types_bar_chart):
                    if is_first or (i%4 == 0):
                        latex.begin_figure(out)
                        is_first = False
                    
                    i += 1
                    latex.write_graphics(out, orf_types_bar_chart, height=0.15)

                    if i%6 == 0:
                        latex.write_caption(out, caption)
                        latex.end_figure(out)
                        latex.clearpage(out)

                else:
                    msg = "Could not find image: {}".format(orf_types_bar_chart)
                    logger.warning(msg)



            if (i > 0) and (i%6) != 0:
                latex.write_caption(out, caption)
                latex.end_figure(out)
                #latex.clearpage(out)

        if i%6 != 0:
            latex.clearpage(out)

    
        # now, if the config file specifies replicates, create figures for those                
        i = 0
        for replicate_name in replicate_names:
            lengths = None
            offsets = None

            caption = "ORF types: {}".format(replicate_name)

            it = itertools.product(grouped_values, chisq_values, filtered_values)

            is_first = True

            for is_grouped, is_chisq, is_filtered in it:

                if is_chisq:
                    f = None
                    rw = None
                else:
                    f = fraction
                    rw = reweighting_iterations

                
                orf_types_bar_chart = filenames.get_orf_types_bar_chart(
                    config['riboseq_data'], 
                    replicate_name, 
                    length=lengths, 
                    offset=offsets, 
                    is_unique=is_unique, 
                    note=out_note_str,
                    image_type=args.image_type,
                    fraction=f,
                    reweighting_iterations=rw,
                    is_grouped=is_grouped,
                    is_chisq=is_chisq,
                    is_filtered=is_filtered
                )

                
                msg = "Looking for image file: {}".format(orf_types_bar_chart)
                logger.debug(msg)

                if os.path.exists(orf_types_bar_chart):
                    if is_first or (i%4 == 0):
                        latex.begin_figure(out)
                        is_first = False
                    
                    i += 1
                    latex.write_graphics(out, orf_types_bar_chart, height=0.15)

                    if i%6 == 0:
                        latex.write_caption(out, caption)
                        latex.end_figure(out)
                        latex.clearpage(out)

                else:
                    msg = "Could not find image: {}".format(orf_types_bar_chart)
                    logger.warning(msg)


            if (i > 0) and (i%6) != 0:
                latex.write_caption(out, caption)
                latex.end_figure(out)
                #latex.clearpage(out)

        if i%6 != 0:
            latex.clearpage(out)


        ### ORF type length distributions
        title = "Predicted ORF type length distributions"
        latex.section(out, title)

        i = 0
        for sample_name in sample_names:
            
            try:
                lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(config,
                                                                               sample_name,
                                                                               is_unique=is_unique,
                                                                               default_params=metagene_options)
            except FileNotFoundError:
                msg = ("Could not parse out lengths and offsets for sample: {}. "
                    "Skipping".format(sample_name))
                logger.error(msg)
                continue
            
            caption = "ORF type length distributions: {}".format(sample_name)

            is_first = True
            it = itertools.product(grouped_values, chisq_values)
            
            for is_grouped, is_chisq in it:

                if is_chisq:
                    f = None
                    rw = None
                else:
                    f = fraction
                    rw = reweighting_iterations

                orf_length_line_graph = filenames.get_orf_length_distribution_line_graph(
                    config['riboseq_data'], 
                    sample_name,
                    length=lengths,
                    offset=offsets, 
                    is_unique=is_unique,
                    note=out_note_str,
                    image_type=args.image_type,
                    fraction=f,
                    reweighting_iterations=rw,
                    is_grouped=is_grouped,
                    is_chisq=is_chisq
                )

                if os.path.exists(orf_length_line_graph):
            
                    if is_first or (i%4 == 0):
                        latex.begin_figure(out)
                        is_first = False
                    
                    i += 1
                    latex.write_graphics(out, orf_length_line_graph, height=0.15)

                    if i%4 == 0:
                        latex.write_caption(out, caption)
                        latex.end_figure(out)
                        latex.clearpage(out)

                else:
                    msg = "Could not find image: {}".format(orf_length_line_graph)
                    logger.debug(msg)


            if (i > 0) and (i%4) != 0:
                latex.write_caption(out, caption)
                latex.end_figure(out)
                #latex.clearpage(out)

        if i%4 != 0:
            latex.clearpage(out)

        # now, if the config file specifies replicates, create figures for those  
        i = 0
        for replicate_name in replicate_names:
            lengths = None
            offsets = None

            caption = "ORF types: {}".format(replicate_name)

            is_first = True
            it = itertools.product(grouped_values, chisq_values)
            
            for is_grouped, is_chisq in it:

                if is_chisq:
                    f = None
                    rw = None
                else:
                    f = fraction
                    rw = reweighting_iterations

                orf_length_line_graph = filenames.get_orf_length_distribution_line_graph(
                    config['riboseq_data'],
                    replicate_name,
                    length=lengths,
                    offset=offsets, 
                    is_unique=is_unique,
                    note=out_note_str,
                    image_type=args.image_type,
                    fraction=f,
                    reweighting_iterations=rw,
                    is_grouped=is_grouped,
                    is_chisq=is_chisq
                )
                
                if os.path.exists(orf_length_line_graph):
            
                    if is_first or (i%4 == 0):
                        latex.begin_figure(out)
                        is_first = False
                    
                    i += 1
                    latex.write_graphics(out, orf_length_line_graph, height=0.15)

                    if i%4 == 0:
                        latex.write_caption(out, caption)
                        latex.end_figure(out)
                        latex.clearpage(out)

                else:
                    msg = "Could not find image: {}".format(orf_length_line_graph)
                    logger.debug(msg)

            if (i > 0) and (i%4) != 0:
                latex.write_caption(out, caption)
                latex.end_figure(out)
                #latex.clearpage(out)

        if i%4 != 0:
            latex.clearpage(out)

        ### ORF type metagene profiles
        if args.show_orf_periodicity:
            title = "Predicted ORF type metagene profiles"
            latex.section(out, title)
            
            i = 0
            for sample_name in sample_names:
                
                try:
                    lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(config,
                                                                                   sample_name,
                                                                                   is_unique=is_unique,
                                                                                   default_params=metagene_options)
                except FileNotFoundError:
                    msg = ("Could not parse out lengths and offsets for sample: {}. "
                        "Skipping".format(sample_name))
                    logger.error(msg)
                    continue
                
                caption = "ORF type metagene profiles: {}".format(sample_name)

                is_first = True

                for is_chisq in chisq_values:

                    if is_chisq:
                        f = None
                        rw = None
                    else:
                        f = fraction
                        rw = reweighting_iterations

                    orf_type_profile_base = filenames.get_orf_type_profile_base(
                        config['riboseq_data'], 
                        sample_name, 
                        length=lengths, offset=offsets, 
                        is_unique=is_unique, 
                        note=out_note_str,
                        fraction=f, reweighting_iterations=rw,
                        is_chisq=is_chisq)

                    it = itertools.product(ribo_utils.orf_types, strands)

                    for orf_type, strand in it:
                        orf_type_profile = filenames.get_orf_type_profile_image(
                            orf_type_profile_base, 
                            orf_type, 
                            strand, 
                            args.image_type
                        )

                        msg = "Looking for image file: {}".format(orf_type_profile)
                        logger.debug(msg)
                        if os.path.exists(orf_type_profile):

                            if is_first or (i%4 == 0):
                                latex.begin_figure(out)
                                is_first = False

                            i += 1
                            latex.write_graphics(out, orf_type_profile, height=0.23)

                            if i%4 == 0:
                                latex.write_caption(out, caption)
                                latex.end_figure(out)
                                latex.clearpage(out)

                        else:
                            msg = "Could not find image: {}".format(orf_type_profile)
                            logger.warning(msg)

                if (i > 0) and (i%4 != 0):
                    latex.write_caption(out, caption)
                    latex.end_figure(out)
                    #latex.clearpage(out)

            if i%4 != 0:
                latex.clearpage(out)

            i = 0
            for replicate_name in replicate_names:
                lengths = None
                offsets = None
                            
                caption = "ORF type metagene profiles: {}".format(replicate_name)
                is_first = True
                for is_chisq in chisq_values:

                    if is_chisq:
                        f = None
                        rw = None
                    else:
                        f = fraction
                        rw = reweighting_iterations

                    orf_type_profile_base = filenames.get_orf_type_profile_base(
                        config['riboseq_data'], 
                        replicate_name, 
                        length=lengths, offset=offsets, 
                        is_unique=is_unique, 
                        note=out_note_str, 
                        fraction=f, reweighting_iterations=rw,
                        is_chisq=is_chisq
                    )

                    it = itertools.product(ribo_utils.orf_types, strands)

                    for orf_type, strand in it:

                        orf_type_profile = filenames.get_orf_type_profile_image(
                            orf_type_profile_base, 
                            orf_type, 
                            strand, 
                            args.image_type
                        )

                        if os.path.exists(orf_type_profile):

                            if is_first or (i%4 == 0):
                                latex.begin_figure(out)
                                is_first = False

                            i += 1
                            latex.write_graphics(out, orf_type_profile, height=0.23)

                            if i % 4 == 0:
                                latex.write_caption(out, caption)
                                latex.end_figure(out)
                                latex.clearpage(out)
                        else:
                            msg = "Could not find image: {}".format(orf_type_profile)
                            logger.debug(msg)
                
                if (i > 0) and (i%4 != 0):
                    latex.write_caption(out, caption)
                    latex.end_figure(out)
                    #latex.clearpage(out)

            if i%4 != 0:
                latex.clearpage(out)

        latex.end_document(out)

    tex_filename = os.path.basename(tex_file)
    latex.compile(args.out, tex_filename)
    
if __name__ == '__main__':
    main()
