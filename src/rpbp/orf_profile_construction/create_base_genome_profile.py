#! /usr/bin/env python3

"""Preprocessing (flexbar, bowtie2, STAR),
create base genome profile.

"""

import os
import sys
import yaml
import logging
import argparse

from pathlib import Path

import pbiotools.utils.bam_utils as bam_utils
import pbiotools.utils.fastx_utils as fastx_utils
import pbiotools.utils.pgrm_utils as pgrm_utils
import pbiotools.misc.logging_utils as logging_utils
import pbiotools.misc.shell_utils as shell_utils
import pbiotools.misc.utils as utils

import rpbp.ribo_utils.filenames as filenames

from rpbp.defaults import (
    default_num_cpus,
    default_mem,
    star_executable,
    star_options,
    flexbar_options,
)

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Creates base genome profile.",
    )

    parser.add_argument("raw_data", help="The raw data file (fastq[.gz])")

    parser.add_argument("config", help="The (yaml) configuration file")

    parser.add_argument(
        "name", help="The name for the dataset, used in the created files"
    )

    parser.add_argument(
        "-p",
        "--num-cpus",
        help="The number of processors to use",
        type=int,
        default=default_num_cpus,
    )

    parser.add_argument(
        "--mem", help="The amount of RAM to request", default=default_mem
    )

    parser.add_argument(
        "-t",
        "--tmp",
        help="""The location for temporary files. If not
        specified, program-specific temp locations are used.""",
        default=None,
    )

    parser.add_argument("--do-not-call", action="store_true")

    parser.add_argument(
        "--overwrite",
        help="""If this flag is present, existing files
        will be overwritten.""",
        action="store_true",
    )

    parser.add_argument(
        "-k",
        "--keep-intermediate-files",
        help="""If this flag is given,
        then all intermediate files will be kept; otherwise, they will be
        deleted. This feature is implemented piecemeal. If the --do-not-call flag
        is given, then nothing will be deleted.""",
        action="store_true",
    )

    logging_utils.add_logging_options(parser)
    pgrm_utils.add_star_options(parser, star_executable)
    pgrm_utils.add_flexbar_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "[create-base-genome-profile]: {}".format(" ".join(sys.argv))
    logger.info(msg)

    config = yaml.load(open(args.config), Loader=yaml.FullLoader)

    # check that all of the necessary programs are callable
    programs = [
        "flexbar",
        args.star_executable,
        "samtools",
        "bowtie2",
        "remove-multimapping-reads",
    ]
    shell_utils.check_programs_exist(programs)

    required_keys = [
        "riboseq_data",
        "ribosomal_index",
        "gtf",
        "genome_base_path",
        "genome_name",
    ]
    utils.check_keys_exist(config, required_keys)

    note = config.get("note", None)
    call = not args.do_not_call
    keep_delete_files = args.keep_intermediate_files or args.do_not_call

    # Step 0: Running flexbar to remove adapter sequences

    raw_data = args.raw_data
    flexbar_target = filenames.get_without_adapters_base(
        config["riboseq_data"], args.name, note=note
    )
    without_adapters = filenames.get_without_adapters_fastq(
        config["riboseq_data"], args.name, note=note
    )

    adapter_seq_str = utils.get_config_argument(
        config, "adapter_sequence", "adapter-seq"
    )
    adapter_file_str = utils.get_config_argument(config, "adapter_file", "adapters")

    # get all options, command line options override defaults
    flexbar_option_str = pgrm_utils.get_final_args(
        flexbar_options, args.flexbar_options
    )

    cmd = "flexbar -r {} -t {} {} {} {} -n {}".format(
        raw_data,
        flexbar_target,
        adapter_seq_str,
        adapter_file_str,
        flexbar_option_str,
        args.num_cpus,
    )
    in_files = [raw_data]
    out_files = [without_adapters]
    file_checkers = {without_adapters: fastx_utils.check_fastq_file}
    shell_utils.call_if_not_exists(
        cmd,
        out_files,
        in_files=in_files,
        file_checkers=file_checkers,
        overwrite=args.overwrite,
        call=call,
    )

    # Step 1: Running bowtie2 to remove rRNA alignments

    out = utils.abspath("dev", "null")  # we do not care about the alignments
    without_rrna = filenames.get_without_rrna_fastq(
        config["riboseq_data"], args.name, note=note
    )
    with_rrna = filenames.get_with_rrna_fastq(
        config["riboseq_data"], args.name, note=note
    )

    cmd = "bowtie2 -p {} --very-fast -x {} -U {} -S {} --un-gz {} --al-gz {}".format(
        args.num_cpus,
        config["ribosomal_index"],
        without_adapters,
        out,
        without_rrna,
        with_rrna,
    )

    in_files = [without_adapters]
    in_files.extend(pgrm_utils.get_bowtie2_index_files(config["ribosomal_index"]))
    out_files = [without_rrna, with_rrna]
    to_delete = [without_adapters]
    file_checkers = {without_rrna: fastx_utils.check_fastq_file}
    shell_utils.call_if_not_exists(
        cmd,
        out_files,
        in_files=in_files,
        file_checkers=file_checkers,
        overwrite=args.overwrite,
        call=call,
        keep_delete_files=keep_delete_files,
        to_delete=to_delete,
    )

    # Step 2: Running STAR to align rRNA-depleted reads to genome

    # STAR standard output
    star_output_prefix = filenames.get_riboseq_bam_base(
        config["riboseq_data"], args.name, note=note
    )
    genome_star_bam = "{}{}".format(star_output_prefix, "Aligned.sortedByCoord.out.bam")
    # Rp-Bp file name
    genome_sorted_bam = filenames.get_riboseq_bam(
        config["riboseq_data"], args.name, note=note
    )

    # get all options, command line options override defaults

    mem_bytes = utils.human2bytes(args.mem)
    star_options["limitBAMsortRAM"] = mem_bytes

    if args.tmp is not None:
        star_tmp_name = str(args.name + "_STARtmp")
        star_tmp_dir = pgrm_utils.create_star_tmp(args.tmp, star_tmp_name)
        star_options["outTmpDir"] = star_tmp_dir

    star_option_str = pgrm_utils.get_final_args(star_options, args.star_options)

    gtf_file = filenames.get_gtf(config)

    cmd = (
        "{} --runThreadN {} --genomeDir {} --sjdbGTFfile {} --readFilesIn {} "
        "{} --outFileNamePrefix {}".format(
            args.star_executable,
            args.num_cpus,
            config["star_index"],
            gtf_file,
            without_rrna,
            star_option_str,
            star_output_prefix,
        )
    )
    in_files = [without_rrna]
    in_files.extend(pgrm_utils.get_star_index_files(config["star_index"]))
    to_delete = [without_rrna]
    # run if any of the two doesn't exist, but only validate genome_star_bam
    # as genome_sorted_bam normally doesn't exist yet
    out_files = [genome_star_bam, genome_sorted_bam]
    file_checkers = {genome_star_bam: bam_utils.check_bam_file}
    shell_utils.call_if_not_exists(
        cmd,
        out_files,
        in_files=in_files,
        file_checkers=file_checkers,
        overwrite=args.overwrite,
        call=call,
        keep_delete_files=keep_delete_files,
        to_delete=to_delete,
    )

    # rename STAR output to that expected by the pipeline
    genome_star_bam = Path(genome_star_bam)
    genome_star_bam.replace(genome_sorted_bam)

    # create the bamtools index
    cmd = "samtools index -b {}".format(genome_sorted_bam)
    shell_utils.check_call(cmd, call=call)

    # check if we want to keep multimappers
    if config.get("keep_riboseq_multimappers", False):
        return

    # remove multimapping reads from the genome file
    tmp_str = ""
    if args.tmp is not None:
        tmp_str = "--tmp {}".format(args.tmp)

    unique_genome_filename = filenames.get_riboseq_bam(
        config["riboseq_data"], args.name, is_unique=True, note=note
    )

    cmd = "remove-multimapping-reads {} {} {}".format(
        genome_sorted_bam, unique_genome_filename, tmp_str
    )

    in_files = [genome_sorted_bam]
    out_files = [unique_genome_filename]
    to_delete = [genome_sorted_bam]
    file_checkers = {unique_genome_filename: bam_utils.check_bam_file}
    shell_utils.call_if_not_exists(
        cmd,
        out_files,
        in_files=in_files,
        file_checkers=file_checkers,
        overwrite=args.overwrite,
        call=call,
        keep_delete_files=keep_delete_files,
        to_delete=to_delete,
    )


if __name__ == "__main__":
    main()
