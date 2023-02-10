#! /usr/bin/env python3

"""This is the main script used to create the reference
genome indices, identify and label the ORFs.

Calls:
    extract-orf-coordinates
    label-orfs
"""

import os
import sys
import yaml
import logging
import argparse

from pathlib import Path

import pbiotools.misc.logging_utils as logging_utils
import pbiotools.misc.shell_utils as shell_utils
import pbiotools.misc.slurm as slurm
import pbiotools.misc.utils as utils

import pbiotools.utils.bed_utils as bed_utils
import pbiotools.utils.pgrm_utils as pgrm_utils

import rpbp.ribo_utils.filenames as filenames

from rpbp.defaults import (
    default_num_cpus,
    default_mem,
    star_executable,
    star_options_genome,
    default_start_codons,
    default_stop_codons,
)

logger = logging.getLogger(__name__)


def get_orfs(gtf, args, config, is_annotated=False, is_de_novo=False):
    """Process a GTF file into its ORFs."""

    call = not args.do_not_call
    chr_name_file = os.path.join(config["star_index"], "chrName.txt")
    chr_name_str = "--chr-name-file {}".format(chr_name_file)

    logging_str = logging_utils.get_logging_options_string(args)
    cpus_str = "--num-cpus {}".format(args.num_cpus)

    # extract a BED12 of the annotated ORFs
    transcript_bed = filenames.get_bed(
        config["genome_base_path"],
        config["genome_name"],
        is_annotated=is_annotated,
        is_de_novo=is_de_novo,
    )

    cmd = "gtf-to-bed12 {} {} {} {} {}".format(
        gtf, transcript_bed, chr_name_str, cpus_str, logging_str
    )
    in_files = [gtf]
    out_files = [transcript_bed]
    shell_utils.call_if_not_exists(
        cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call
    )

    # extract the transcript fasta
    transcript_fasta = filenames.get_transcript_fasta(
        config["genome_base_path"],
        config["genome_name"],
        is_annotated=is_annotated,
        is_de_novo=is_de_novo,
    )

    cmd = "extract-bed-sequences {} {} {} {}".format(
        transcript_bed, config["fasta"], transcript_fasta, logging_str
    )
    in_files = [transcript_bed, config["fasta"]]
    out_files = [transcript_fasta]
    shell_utils.call_if_not_exists(
        cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call
    )

    # extract ORFs from the transcripts using genomic coordinates
    orfs_genomic = filenames.get_orfs(
        config["genome_base_path"],
        config["genome_name"],
        note=config.get("orf_note"),
        is_annotated=is_annotated,
        is_de_novo=is_de_novo,
    )

    start_codons_str = utils.get_config_argument(
        config, "start_codons", default=default_start_codons
    )

    stop_codons_str = utils.get_config_argument(
        config, "stop_codons", default=default_stop_codons
    )

    cmd = "extract-orf-coordinates {} {} {} {} {} {} {}".format(
        transcript_bed,
        transcript_fasta,
        orfs_genomic,
        cpus_str,
        start_codons_str,
        stop_codons_str,
        logging_str,
    )
    in_files = [transcript_fasta, transcript_bed]
    out_files = [orfs_genomic]
    shell_utils.call_if_not_exists(
        cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call
    )

    # write the ORF exons, used to label the ORFs
    exons_file = filenames.get_exons(
        config["genome_base_path"],
        config["genome_name"],
        note=config.get("orf_note"),
        is_annotated=is_annotated,
        is_de_novo=is_de_novo,
    )

    cmd = "split-bed12-blocks {} {} --num-cpus {} {}".format(
        orfs_genomic, exons_file, args.num_cpus, logging_str
    )
    in_files = [orfs_genomic]
    out_files = [exons_file]
    shell_utils.call_if_not_exists(
        cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call
    )

    # label the ORFs
    labeled_orfs = filenames.get_labels(
        config["genome_base_path"],
        config["genome_name"],
        note=config.get("orf_note"),
        is_annotated=is_annotated,
        is_de_novo=is_de_novo,
    )

    annotated_bed = filenames.get_bed(
        config["genome_base_path"],
        config["genome_name"],
        is_annotated=True,
    )

    orf_exons_str = "--orf-exons {}".format(exons_file)

    de_novo_str = ""
    if is_de_novo:
        de_novo_str = '--label-prefix "novel_" --filter --nonoverlapping-label "novel"'

    cmd = "label-orfs {} {} {} {} {} {} {}".format(
        annotated_bed,
        orfs_genomic,
        labeled_orfs,
        orf_exons_str,
        de_novo_str,
        logging_str,
        cpus_str,
    )
    in_files = [annotated_bed, orfs_genomic, exons_file]
    #  ** this function overwrites the input file `orfs_genomic`
    out_files = [labeled_orfs]
    shell_utils.call_if_not_exists(
        cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call
    )


def get_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Prepare genome indices and annotations "
        "to use with the Rp-Bp profile construction/periodicity estimation "
        "and ORF discovery pipeline.",
    )

    parser.add_argument("config", help="A YAML configuration file.")

    parser.add_argument(
        "--overwrite",
        help="Overwrite existing output.",
        action="store_true",
    )

    pgrm_utils.add_star_options(parser, star_executable)
    slurm.add_sbatch_options(parser, num_cpus=default_num_cpus, mem=default_mem)
    logging_utils.add_logging_options(parser)

    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "[prepare-rpbp-genome]: {}".format(" ".join(sys.argv))
    logger.info(msg)

    config = yaml.load(open(args.config), Loader=yaml.FullLoader)

    # check required callable programs, config keys and files
    programs = [
        "extract-orf-coordinates",
        "label-orfs",
        "bowtie2-build-s",
        "split-bed12-blocks",
        "gtf-to-bed12",
        args.star_executable,
    ]
    shell_utils.check_programs_exist(programs)

    required_keys = [
        "genome_base_path",
        "genome_name",
        "gtf",
        "fasta",
        "ribosomal_fasta",
        "ribosomal_index",
        "star_index",
    ]
    utils.check_keys_exist(config, required_keys)

    files = [config["gtf"], config["fasta"], config["ribosomal_fasta"]]
    if "de_novo_gtf" in config:
        files += [config["de_novo_gtf"]]
    utils.check_files_exist(files, source="prepare-rpbp-genome")

    if args.use_slurm:
        cmd = " ".join(sys.argv)
        slurm.check_sbatch(cmd, args=args)
        return

    call = not args.do_not_call

    # the rRNA index
    cmd = "bowtie2-build-s {} {}".format(
        config["ribosomal_fasta"], config["ribosomal_index"]
    )

    in_files = [config["ribosomal_fasta"]]
    out_files = pgrm_utils.get_bowtie2_index_files(config["ribosomal_index"])
    shell_utils.call_if_not_exists(
        cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call
    )

    # the STAR index
    mem = utils.human2bytes(args.mem)
    star_option_str = pgrm_utils.get_final_args(star_options_genome, args.star_options)
    cmd = (
        "{} --runMode genomeGenerate --genomeDir {} --genomeFastaFiles {} "
        "--runThreadN {} --limitGenomeGenerateRAM {} {}".format(
            args.star_executable,
            config["star_index"],
            config["fasta"],
            args.num_cpus,
            mem,
            star_option_str,
        )
    )

    in_files = [config["fasta"]]
    out_files = pgrm_utils.get_star_index_files(config["star_index"])
    shell_utils.call_if_not_exists(
        cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call
    )

    # get the ORFs
    get_orfs(config["gtf"], args, config, is_annotated=True, is_de_novo=False)

    # we will use these files later in the pipeline
    annotated_orfs = filenames.get_orfs(
        config["genome_base_path"],
        config["genome_name"],
        note=config.get("orf_note"),
        is_annotated=True,
        is_de_novo=False,
    )

    orfs_genomic = filenames.get_orfs(
        config["genome_base_path"], config["genome_name"], note=config.get("orf_note")
    )

    annotated_exons_file = filenames.get_exons(
        config["genome_base_path"],
        config["genome_name"],
        note=config.get("orf_note"),
        is_annotated=True,
        is_de_novo=False,
    )

    exons_file = filenames.get_exons(
        config["genome_base_path"], config["genome_name"], note=config.get("orf_note")
    )

    annotated_labeled_orfs = filenames.get_labels(
        config["genome_base_path"],
        config["genome_name"],
        note=config.get("orf_note"),
        is_annotated=True,
        is_de_novo=False,
    )

    labeled_orfs = filenames.get_labels(
        config["genome_base_path"], config["genome_name"], note=config.get("orf_note")
    )

    # now, check if we have a de novo assembly
    if "de_novo_gtf" in config:
        get_orfs(
            config["de_novo_gtf"], args, config, is_annotated=False, is_de_novo=True
        )

        # we need to concat the ORF and exon files
        de_novo_orfs = filenames.get_orfs(
            config["genome_base_path"],
            config["genome_name"],
            note=config.get("orf_note"),
            is_annotated=False,
            is_de_novo=True,
        )

        orfs_files = [annotated_orfs, de_novo_orfs]

        orfs_files_str = " ".join(orfs_files)
        msg = "Concatenating files. Output file: {}; Input files: {}".format(
            orfs_genomic, orfs_files_str
        )
        logger.info(msg)

        if call:
            concatenated_bed = bed_utils.concatenate(orfs_files, sort_bed=True)
            concatenated_bed["orf_num"] = range(len(concatenated_bed))
            additional_columns = ["orf_num", "orf_len"]
            fields = bed_utils.bed12_field_names + additional_columns
            bed_utils.write_bed(concatenated_bed[fields], orfs_genomic)
        else:
            msg = "Skipping concatenation due to --call value"
            logger.info(msg)

        de_novo_exons_file = filenames.get_exons(
            config["genome_base_path"],
            config["genome_name"],
            note=config.get("orf_note"),
            is_annotated=False,
            is_de_novo=True,
        )

        exons_files = [annotated_exons_file, de_novo_exons_file]

        exons_files_str = " ".join(exons_files)
        msg = "Concatenating files. Output file: {}; Input files: {}".format(
            exons_file, exons_files_str
        )
        logger.info(msg)

        if call:
            concatenated_bed = bed_utils.concatenate(exons_files, sort_bed=True)
            fields = bed_utils.bed6_field_names + ["exon_index", "transcript_start"]
            bed_utils.write_bed(concatenated_bed[fields], exons_file)
        else:
            msg = "Skipping concatenation due to --call value"
            logger.info(msg)

        de_novo_labeled_orfs = filenames.get_labels(
            config["genome_base_path"],
            config["genome_name"],
            note=config.get("orf_note"),
            is_annotated=False,
            is_de_novo=True,
        )

        label_files = [annotated_labeled_orfs, de_novo_labeled_orfs]

        label_files_str = " ".join(label_files)
        msg = "Concatenating files. Output file: {}; Input files: {}".format(
            labeled_orfs, label_files_str
        )
        logger.info(msg)

        if call:
            # not a BED file
            concatenated_bed = bed_utils.concatenate(label_files, sort_bed=False)
            bed_utils.write_bed(concatenated_bed, labeled_orfs)
        else:
            msg = "Skipping concatenation due to --call value"
            logger.info(msg)

        # concat the annotations to inform STAR
        gtf_file = filenames.get_gtf(config)
        star_files_str = " ".join([config["gtf"], config["de_novo_gtf"]])
        msg = "Concatenating files. Output file: {}; Input files: {}".format(
            gtf_file, star_files_str
        )
        logger.info(msg)

        cmd = "awk '!/^#/' {} {} > {}".format(
            config["gtf"], config["de_novo_gtf"], gtf_file
        )
        in_files = [config["gtf"], config["de_novo_gtf"]]
        out_files = [gtf_file]
        shell_utils.call_if_not_exists(
            cmd, out_files, in_files=in_files, overwrite=args.overwrite, call=call
        )

    else:
        # if we do not have a de novo assembly, rename the files
        annotated_orfs = Path(annotated_orfs)
        annotated_orfs.replace(orfs_genomic)

        annotated_exons_file = Path(annotated_exons_file)
        annotated_exons_file.replace(exons_file)

        annotated_labeled_orfs = Path(annotated_labeled_orfs)
        annotated_labeled_orfs.replace(labeled_orfs)


if __name__ == "__main__":
    main()
