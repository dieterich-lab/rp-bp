#! /usr/bin/env python3

import argparse
import logging
import os
import yaml

import pbiotools.misc.logging_utils as logging_utils
import pbiotools.misc.shell_utils as shell_utils
import pbiotools.misc.utils as utils
import pbiotools.misc.slurm as slurm

import rpbp.ribo_utils.utils as ribo_utils
import rpbp.ribo_utils.filenames as filenames

logger = logging.getLogger(__name__)

default_num_procs = 2
default_note = None

default_peptide_separator = "\t"
default_peptide_filter_field = "PEP"
default_peptide_filter_value = 0.1


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script identifies the orf peptide matches for all samples in "
        "a project.",
    )
    parser.add_argument("config", help="The (yaml) config file")

    parser.add_argument(
        "--peptide-filter-field",
        help="The field to use for " "filtering the peptides from MaxQuant",
        default=default_peptide_filter_field,
    )

    parser.add_argument(
        "--peptide-filter-value",
        help="All peptides with a value "
        "greater than the filter value will be removed",
        type=float,
        default=default_peptide_filter_value,
    )

    parser.add_argument(
        "--peptide-separator",
        help="The separator in the " "peptide file",
        default=default_peptide_separator,
    )

    parser.add_argument(
        "--note",
        help="If this option is given, it will be used in "
        "the output filenames.\n\nN.B. This REPLACES the note in the config file.",
        default=default_note,
    )

    slurm.add_sbatch_options(parser)
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    logging_str = logging_utils.get_logging_options_string(args)

    config = yaml.load(open(args.config), Loader=yaml.FullLoader)

    programs = ["get-orf-peptide-matches"]
    shell_utils.check_programs_exist(programs)

    required_keys = [
        "peptide_files",
        "peptide_cell_type_analysis",
        "riboseq_data",
        "riboseq_samples",
    ]
    utils.check_keys_exist(config, required_keys)

    note_str = config.get("note", None)
    out_note_str = note_str

    if args.note is not None and len(args.note) > 0:
        out_note_str = args.note

    args_dict = vars(args)

    peptide_filter_field_str = utils.get_config_argument(
        args_dict, "peptides_filter_field"
    )
    peptide_filter_value_str = utils.get_config_argument(
        args_dict, "peptides_filter_value"
    )
    peptide_separator_str = utils.get_config_argument(args_dict, "peptide_separator")

    num_cpus_str = utils.get_config_argument(args_dict, "num_cpus")

    cell_types = ribo_utils.get_riboseq_cell_type_samples(config)
    for cell_type, peptide_files in config["peptide_cell_type_analysis"].items():
        if cell_type not in cell_types:
            msg = (
                "Could not find cell_type specification. Please check the config "
                "file: {}".format(cell_type)
            )
            logger.warning(msg)
            continue

        cell_type_protein = filenames.get_riboseq_cell_type_protein(
            config["riboseq_data"], cell_type, is_filtered=True, note=note_str
        )

        if not os.path.exists(cell_type_protein):
            msg = "Could not find cell_type protein fasta. Skipping: {}".format(
                cell_type_protein
            )
            logger.warning(msg)
            continue

        for peptide_file in peptide_files:
            if peptide_file not in config["peptide_files"]:
                msg = (
                    "Could not find peptide_file specification. Please check "
                    "the config file: {}".format(peptide_file)
                )
                logger.warning(msg)
                continue

            peptide_txt_file = config["peptide_files"][peptide_file]

            if not os.path.exists(peptide_txt_file):
                msg = "Could not find peptide.txt file. Skipping: {}".format(
                    peptide_txt_file
                )
                logger.warning(msg)
                continue

            peptide_matches = filenames.get_riboseq_peptide_matches(
                config["riboseq_data"],
                cell_type,
                peptide_file,
                is_filtered=True,
                note=out_note_str,
            )

            cmd = "get-orf-peptide-matches {} {} {} {} {} {} {} {}".format(
                cell_type_protein,
                peptide_txt_file,
                peptide_matches,
                num_cpus_str,
                peptide_filter_field_str,
                peptide_filter_value_str,
                peptide_separator_str,
                logging_str,
            )

            slurm.check_sbatch(cmd, args=args)


if __name__ == "__main__":
    main()
