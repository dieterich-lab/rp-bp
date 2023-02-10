#! /usr/bin/env python3

import argparse
import logging
import shlex
import sys
import re
import yaml
import json
import tempfile
import datetime

from pathlib import Path
from collections import defaultdict

import pandas as pd
import numpy as np

import pbiotools.misc.logging_utils as logging_utils
import pbiotools.misc.parallel as parallel
import pbiotools.misc.shell_utils as shell_utils
import pbiotools.utils.bed_utils as bed_utils
import pbiotools.misc.utils as utils

import rpbp.ribo_utils.utils as ribo_utils
import rpbp.ribo_utils.filenames as filenames

from rpbp.defaults import (
    default_num_cpus,
    metagene_options,
    orf_type_colors,
    orf_type_name_map,
    orf_type_labels,
)

logger = logging.getLogger(__name__)

# Supplement in xlsx format (Phase I Table S2, sheet 2)
# and single-study Ribo-seq ORFs (Table S3, sheet 3).
standardized_orfs_url = (
    "https://static-content.springer.com/"
    "esm/art%3A10.1038%2Fs41587-022-01369-0/MediaObjects/"
    "41587_2022_1369_MOESM2_ESM.xlsx"
)

# appris_url = "https://apprisws.bioinfo.cnio.es/pub/current_release/datafiles"
# appris_file = "appris_data.appris.txt"
# appris_organisms = ["bos_taurus", "caenorhabditis_elegans", "danio_rerio",
# "drosophila_melanogaster", "gallus_gallus", "homo_sapiens",
# "macaca_mulatta", "mus_musculus", "rattus_norvegicus",
# "sus_scrofa"]

MERGE_FIELDS = [
    "seqname",
    "start",
    "end",
    "strand",
    "num_exons",
    "exon_lengths",
    "exon_genomic_relative_starts",
]

fields = [
    "condition",
    "bayes_factor_mean",
    "bayes_factor_var",
    "x_1_sum",
    "x_2_sum",
    "x_3_sum",
    "orf_num",
    "orf_len",
    "orf_type",
    "biotype",
    "transcripts",
    "gene_id",
    "gene_name",
    "gene_biotype",
]

FIELD_NAMES = bed_utils.bed12_field_names + fields[1:8]
EXT_FIELD_NAMES = bed_utils.bed12_field_names + fields


def get_predictions_file(
    name,
    is_sample,
    note,
    fraction,
    reweighting_iterations,
    is_unique,
    is_filtered,
    config,
    ftype="predictions",
):

    if is_sample:
        try:
            lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(
                config, name, is_unique=is_unique, default_params=metagene_options
            )
        except FileNotFoundError:
            msg = (
                "Could not find periodic lengths and offset file."
                f"Skipping. name: {name}."
            )
            logger.error(msg)
            return None
    else:
        lengths, offsets = None, None

    if ftype == "predictions":
        filen = filenames.get_riboseq_predicted_orfs(
            config["riboseq_data"],
            name,
            length=lengths,
            offset=offsets,
            is_unique=is_unique,
            note=note,
            fraction=fraction,
            reweighting_iterations=reweighting_iterations,
            is_filtered=is_filtered,
        )
    elif ftype == "profiles":
        filen = filenames.get_riboseq_profiles(
            config["riboseq_data"],
            name,
            length=lengths,
            offset=offsets,
            is_unique=is_unique,
            note=note,
            fraction=fraction,
            reweighting_iterations=reweighting_iterations,
        )
    elif ftype == "base":
        filen = filenames.get_orf_type_profile_base(
            config["riboseq_data"],
            name,
            length=lengths,
            offset=offsets,
            is_unique=is_unique,
            note=note,
            fraction=fraction,
            reweighting_iterations=reweighting_iterations,
        )
    else:
        msg = f"File type {ftype} unrecognized!"
        logger.critical(msg)

    if not Path(filen).is_file() and not ftype == "base":
        msg = f"Could not find predicted ORFs. name: {name}, file: {filen}."
        raise FileNotFoundError(msg)

    return filen


def chrom_match(match, string):
    return bool(re.match(match, string))


def cut_it(x, seqname, karyotype, args):
    size = int(karyotype[seqname])
    bins = np.arange(0, size + args.circos_bin_width, args.circos_bin_width)
    if (bins[-1] - size) < args.circos_bin_width / 2:
        bins[-1] = size
    else:
        bins[-2] = size
        bins = np.delete(bins, -1)
    return pd.cut(x["start"].values, bins, right=False)


def get_circos_graph(orfs, sub_folder, config, args):

    orf_types = orfs.orf_type.unique()
    df = orfs.drop_duplicates(subset=bed_utils.bed12_field_names).copy()
    df.rename(columns={"seqname": "block_id"}, inplace=True)
    df.reset_index(drop=True, inplace=True)

    # filter chromosomes to show
    match = r"|".join(args.circos_show_chroms)
    m_block = df["block_id"].apply(lambda x: chrom_match(match, x))
    df = df[m_block]

    filen = Path(config["star_index"], "chrNameLength.txt")
    with open(filen) as f:
        karyotype = dict(x.rstrip().split(None, 1) for x in f)
    karyotype = {k: v for k, v in karyotype.items() if k in df.block_id.unique()}

    vals = (
        df.groupby("block_id")
        .apply(lambda x: cut_it(x, x.name, karyotype, args))
        .values
    )
    df["range"] = utils.flatten_lists(vals)
    df["start"] = [v.left for v in df["range"].values]
    df["end"] = [
        v.right - 1 for v in df["range"].values
    ]  # we're one off for the last bin...
    fields = ["block_id", "orf_type", "start", "end"]
    records = (
        df.groupby(fields)["id"]
        .count()
        .reset_index(name="value")
        .to_dict(orient="records")
    )
    circos_graph_data = defaultdict(list)
    for orf_type in orf_types:
        for record in records:
            if not record["orf_type"] == orf_type:
                continue
            clean = {k: v for k, v in record.items() if not k == "orf_type"}
            circos_graph_data[f"histogram_{orf_type}"].append(clean)

    colors = ["#ededed", "#333333"]
    layout = []
    for idx, (seqname, size) in enumerate(karyotype.items()):
        label = seqname if "chr" in seqname else f"chr{seqname}"
        record = {
            "id": seqname,
            "label": label,
            "color": colors[idx % 2],
            "len": int(size),
        }
        layout.append(record)
    circos_graph_data["genome"] = layout

    filen = f"{config['genome_name']}.circos_graph_data.json"
    filen = Path(config["riboseq_data"], sub_folder, filen)
    with open(filen, "w") as f:
        json.dump(circos_graph_data, f)

    return filen.as_posix()


def add_data(
    name,
    sample_name_map,
    is_sample,
    note,
    fraction,
    reweighting_iterations,
    is_unique,
    is_filtered,
    config,
):

    orfs_file = get_predictions_file(
        name,
        is_sample,
        note,
        fraction,
        reweighting_iterations,
        is_unique,
        is_filtered,
        config,
    )
    orfs = bed_utils.read_bed(orfs_file, low_memory=False)
    orfs = orfs[FIELD_NAMES]
    orfs["condition"] = sample_name_map[name]

    return orfs


def get_bed_blocks(row):

    # convert genomic coordinates to bed blocks
    # ad-hoc for standardized ORFs

    starts = row.starts.split(";")
    ends = row.ends.split(";")
    # coordinates include stop codon
    # a few ORFs have a stop codon overlapping 2 exons
    # since we are getting rid of it, we can discard these blocks
    if row.strand == "+":
        if abs(int(ends[-1]) - int(starts[-1])) <= 3:
            starts = starts[:-1]
            ends = ends[:-1]
        else:
            stop = int(ends[-1])
            stop -= 3
            ends[-1] = str(stop)
    else:
        if abs(int(ends[0]) - int(starts[0])) <= 3:
            starts = starts[1:]
            ends = ends[1:]
        else:
            stop = int(starts[0])
            stop += 3
            starts[0] = str(stop)

    coords = [f"{start}-{end}" for start, end in zip(starts, ends)]
    coords = ",".join(coords)
    coords = bed_utils.convert_genomic_coords_to_bed_blocks(coords)
    # we now subtract 1 from the start because BED is base-0
    coords["start"] = int(starts[0]) - 1
    coords["end"] = int(ends[-1])

    return coords


def get_standardized_orfs(filen, sheet):

    # ad-hoc for standardized ORFs

    colmap = [
        {"chrm": "seqname", "orf_name": "PHASE I ORFs"},
        {"chrm": "seqname", "orf_name": "SS ORFs"},
    ]
    col = colmap[sheet - 2]["orf_name"]
    standardized_fields = MERGE_FIELDS + [col]

    standardized_orfs = pd.read_excel(filen, sheet_name=sheet)
    blocks = standardized_orfs.apply(get_bed_blocks, axis=1)
    blocks = pd.DataFrame(blocks.to_list())
    standardized_orfs = standardized_orfs.join(blocks)
    standardized_orfs.rename(columns=colmap[sheet - 2], inplace=True)

    return col, standardized_orfs[standardized_fields]


def _create_figures(name_pretty_name_is_sample, config, args):

    name, pretty_name, is_sample = name_pretty_name_is_sample

    is_unique = not config.get("keep_riboseq_multimappers", False)
    note = config.get("note", None)
    fraction = config.get("smoothing_fraction", None)
    reweighting_iterations = config.get("smoothing_reweighting_iterations", None)

    is_filtered = not args.use_unfiltered

    logging_str = logging_utils.get_logging_options_string(args)
    image_type_str = "--image-type {}".format(args.image_type)

    orfs = get_predictions_file(
        name,
        is_sample,
        note,
        fraction,
        reweighting_iterations,
        is_unique,
        is_filtered,
        config,
    )

    profiles = get_predictions_file(
        name,
        is_sample,
        note,
        fraction,
        reweighting_iterations,
        is_unique,
        is_filtered,
        config,
        ftype="profiles",
    )

    orf_type_profile_base = get_predictions_file(
        name,
        is_sample,
        note,
        fraction,
        reweighting_iterations,
        is_unique,
        is_filtered,
        config,
        ftype="base",
    )

    title_str = "{}, Rp-Bp".format(pretty_name)
    title_str = shlex.quote(title_str)
    title_str = "--title {}".format(title_str)

    strand = "+"
    orf_type_profiles_forward = [
        filenames.get_orf_type_profile_image(
            orf_type_profile_base, orf_type, strand, args.image_type
        )
        for orf_type in orf_type_name_map.values()
    ]

    strand = "-"
    orf_type_profiles_reverse = [
        filenames.get_orf_type_profile_image(
            orf_type_profile_base, orf_type, strand, args.image_type
        )
        for orf_type in orf_type_name_map.values()
    ]

    cmd = "visualize-orf-type-metagene-profiles {} {} {} {} {} {} {}".format(
        args.config,
        orfs,
        profiles,
        orf_type_profile_base,
        title_str,
        image_type_str,
        logging_str,
    )

    in_files = [orfs]
    out_files = orf_type_profiles_forward + orf_type_profiles_reverse
    shell_utils.call_if_not_exists(
        cmd, out_files, in_files=in_files, overwrite=args.overwrite
    )


def create_all_figures(config, sample_name_map, condition_name_map, args):

    is_sample = True
    sample_names = sorted(config["riboseq_samples"].keys())
    samples = [(name, sample_name_map[name], is_sample) for name in sample_names]

    conditions = []
    if not args.no_replicates:
        is_sample = False
        replicate_names = sorted(ribo_utils.get_riboseq_replicates(config).keys())
        conditions = [
            (name, condition_name_map[name], is_sample) for name in replicate_names
        ]

    all_names = samples + conditions
    parallel.apply_parallel_iter(
        all_names, default_num_cpus, _create_figures, config, args, progress_bar=True
    )


def get_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Summarizes the ORF prediction step "
        "and prepare data for the web application.",
    )

    parser.add_argument(
        "config",
        help="A YAML configuration file. " "The same used to run the pipeline.",
    )

    # post-hoc filters - other filters such as ORF min length, etc. must be
    # set in the config before running the pipeline
    parser.add_argument(
        "--min-samples",
        help="An ORF is filtered out if not predicted in at "
        "least ``--min-samples`` number of samples. By default "
        "all ORFs are kept. This is ignored if merged replicates "
        "are included in the output.",
        type=int,
        default=1,
    )

    parser.add_argument(
        "-k",
        "--keep-other",
        help="Include ORFs labeled as *other*, if present. "
        "They are discarded by default.",
        action="store_true",
    )

    parser.add_argument(
        "--no-replicates",
        help="If Rp-Bp was run with ``--merge-replicates``, "
        "predictions from merged replicates are included by "
        "default, unless this flag is present.",
        required="--min-samples" in sys.argv,
        action="store_true",
    )

    parser.add_argument(
        "--use-unfiltered",
        help="Use the *unfiltered* ORF predictions. "
        "Unless Rp-Bp was run with ``--write-unfiltered``, "
        "these will not be available. By default, the "
        "*filtered* predictions are used.",
        action="store_true",
    )

    # display
    parser.add_argument(
        "--use-name-maps",
        help="Use ``riboseq_sample_name_map`` and ``riboseq_condition_name_map`` "
        "from the config. Do not use this flag when preparing results for "
        "the dashboard, mapping is done in the app.",
        action="store_true",
    )

    # extra
    parser.add_argument(
        "--match-standardized-orfs",
        help="Add matching ORF names from "
        "https://doi.org/10.1038/s41587-022-01369-0 (Human data)",
        action="store_true",
    )

    parser.add_argument(
        "--circos-bin-width",
        help="Bin width for counting ORF predictions along chromosomes. "
        "Same bin width for all chromosomes. Last bin adjusted ad hoc "
        "to chromosome size.",
        type=int,
        default=10000000,
    )

    parser.add_argument(
        "--circos-show-chroms",
        help="A list of chromosomes for which predictions are shown. "
        "By default, only numbered chromosomes (and X/x, Y/y) are shown, "
        "to avoid cluttering the figure. Use this option for organisms "
        "with a different nomenclature, or to show additional chromosomes "
        "or contigs. ",
        nargs="+",
        default=["\d", "X", "x", "Y", "y"],
    )

    # subparser = parser.add_subparsers()
    # parser_appris = subparser.add_parser('APPRIS',
    # help='Add APPRIS Annotation to transcripts.'
    # )

    # parser_appris.add_argument('appris-organism', help="Organism.",
    # type=str, choices=appris_organisms
    # )

    # parser_appris.add_argument('appris-assembly-version', help="""Assemble version/
    # Gene Dataset. Must be a valid entry from """,
    # type=str)

    # keep some figures from the old reporting tool for debugging
    parser.add_argument(
        "--show-orf-periodicity",
        help="Create bar charts showing the periodicity "
        "of each ORF type (not for the app).",
        action="store_true",
    )

    parser.add_argument(
        "--image-type",
        help="Format for ``--show-orf-periodicity``.",
        default="eps",
    )

    parser.add_argument(
        "--overwrite",
        help="Overwrite existing output.",
        action="store_true",
    )

    logging_utils.add_logging_options(parser)

    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "[summarize-rpbp-predictions]: {}".format(" ".join(sys.argv))
    logger.info(msg)

    config = yaml.load(open(args.config), Loader=yaml.FullLoader)

    programs = [
        "visualize-orf-type-metagene-profiles",
    ]
    shell_utils.check_programs_exist(programs)

    required_keys = [
        "riboseq_data",
        "riboseq_samples",
        "genome_base_path",
        "genome_name",
        "star_index",
    ]
    utils.check_keys_exist(config, required_keys)

    # create directory for summary data
    sub_folder = Path("analysis", "rpbp_predictions")
    Path(config["riboseq_data"], sub_folder).mkdir(parents=True, exist_ok=True)

    # nomenclature
    project = config.get("project_name", "rpbp")
    note = config.get("note", None)
    orf_note = config.get("orf_note", None)
    is_unique = not config.get("keep_riboseq_multimappers", False)
    # defaults are unused in file names
    fraction = config.get("smoothing_fraction", None)
    reweighting_iterations = config.get("smoothing_reweighting_iterations", None)

    is_filtered = not args.use_unfiltered

    ext_field_names = EXT_FIELD_NAMES.copy()

    # option handling
    if not args.no_replicates:
        if "riboseq_biological_replicates" not in config:
            args.no_replicates = True
        else:
            if config["riboseq_biological_replicates"] is None:
                args.no_replicates = True
        if args.no_replicates:
            msg = (
                "Did not find 'riboseq_biological_replicates' key "
                "in config, setting [--no-replicates] to True."
            )
            logger.warning(msg)
        else:  # checking if they were actually run
            riboseq_replicates = ribo_utils.get_riboseq_replicates(config)
            for replicate in riboseq_replicates.keys():
                try:
                    _ = get_predictions_file(
                        replicate,
                        False,
                        note,
                        fraction,
                        reweighting_iterations,
                        is_unique,
                        is_filtered,
                        config,
                    )
                except:
                    args.no_replicates = True
            if args.no_replicates:
                msg = (
                    "The 'riboseq_biological_replicates' key was found, "
                    "but some predictions were missing, setting "
                    "[--no-replicates] to True."
                )
                logger.warning(msg)

    # Collecting all predictions
    msg = "Parsing predictions for samples"
    logger.info(msg)

    cfg = {}
    if args.use_name_maps:
        cfg = config
    sample_name_map = ribo_utils.get_sample_name_map(cfg)
    condition_name_map = ribo_utils.get_condition_name_map(cfg)

    is_sample = True
    all_predictions = parallel.apply_iter_simple(
        config["riboseq_samples"].keys(),
        add_data,
        sample_name_map,
        is_sample,
        note,
        fraction,
        reweighting_iterations,
        is_unique,
        is_filtered,
        config,
    )

    if not args.no_replicates:
        msg = "Parsing predictions for merged replicates"
        logger.info(msg)

        is_sample = False
        condition_predictions = parallel.apply_iter_simple(
            riboseq_replicates,
            add_data,
            condition_name_map,
            is_sample,
            note,
            fraction,
            reweighting_iterations,
            is_unique,
            is_filtered,
            config,
        )
        all_predictions = all_predictions + condition_predictions

    orfs = pd.concat(all_predictions)

    # Adding ORF labels and transcript/gene info
    msg = "Adding information to ORFs"
    logger.info(msg)

    transcript_bed = filenames.get_bed(
        config["genome_base_path"], config["genome_name"], is_annotated=True
    )
    cols = ["id", "biotype", "gene_id", "gene_name", "gene_biotype"]
    bed_df = bed_utils.read_bed(transcript_bed, low_memory=False)[cols]
    bed_df.rename(columns={"id": "transcript_id"}, inplace=True)

    # add info for de novo
    if "de_novo_gtf" in config:
        transcript_bed_dn = filenames.get_bed(
            config["genome_base_path"],
            config["genome_name"],
            is_annotated=False,
            is_de_novo=True,
        )
        cols = ["id", "biotype", "gene_id", "gene_name", "gene_biotype"]
        bed_df_dn = bed_utils.read_bed(transcript_bed_dn, low_memory=False)[cols]
        bed_df_dn.rename(columns={"id": "transcript_id"}, inplace=True)
        bed_df = pd.concat([bed_df, bed_df_dn])

    labeled_orfs = filenames.get_labels(
        config["genome_base_path"], config["genome_name"], note=orf_note
    )
    cols = ["id", "orf_type", "transcripts"]
    labels_df = bed_utils.read_bed(labeled_orfs)[cols]

    orfs = pd.merge(orfs, labels_df, on="id", how="left")
    # pick assigned transcript to annotate ORF
    orfs[["transcript_id", "others"]] = orfs["transcripts"].str.split(
        ",", 1, expand=True
    )
    orfs = pd.merge(orfs, bed_df, on="transcript_id", how="left")

    # Applying post-hoc filtering
    msg = "Applying post-hoc filtering"
    logger.info(msg)

    remaining_labels = ", ".join(orfs.orf_type.unique())
    msg = (
        f"Found these ORF labels: {remaining_labels} "
        f"with frequency\n{orfs.orf_type.value_counts().to_string()}"
    )
    logger.info(msg)

    # add "group" label for processing
    for orf_type, labels in orf_type_labels.items():
        types = [orf_type_name_map[label] for label in labels]
        orfs.loc[orfs["orf_type"].isin(types), "orf_category"] = orf_type
    # just to make sure...
    remove_m = orfs["orf_type"].isna()
    if not args.keep_other:
        other_m = orfs["orf_category"] == "other"
        remove_m = remove_m | other_m
    orfs = orfs[~remove_m]

    remaining_labels = ", ".join(orfs.orf_type.unique())
    msg = (
        f"Remaining ORF labels: {remaining_labels} "
        f"with frequency\n{orfs.orf_type.value_counts().to_string()}"
    )
    logger.info(msg)

    if args.no_replicates and args.min_samples > 1:
        num_orfs = len(orfs["id"].unique())
        orfs = orfs.groupby("id").filter(lambda x: len(x) >= args.min_samples)
        discarded_orfs = num_orfs - len(orfs["id"].unique())
        msg = (
            f"Using [--no-replicates] and [--min-samples {args.min_samples}]: "
            f"Removing {discarded_orfs} ORFs."
        )
        logger.info(msg)

    # Adding standardized ORFs - hard coded
    # based on published data Table S2 and S3
    if args.match_standardized_orfs:
        msg = "Adding standardized ORFs (Human)"
        logger.info(msg)
        local_filename = Path(standardized_orfs_url).name
        with tempfile.TemporaryDirectory() as tmpdirname:
            filen = shell_utils.download_file(
                standardized_orfs_url, local_filename=Path(tmpdirname, local_filename)
            )
            for sheet in range(2, 4):
                field, standardized_orfs = get_standardized_orfs(filen, sheet)
                orfs = pd.merge(orfs, standardized_orfs, on=MERGE_FIELDS, how="left")

                ext_field_names += [field]

    # writing ORFs to disk
    # sort on the chrom field, and then on the chromStart field.
    orfs.sort_values(["seqname", "start"], ascending=[True, True], inplace=True)
    orfs_filen = filenames.get_riboseq_predicted_orfs(
        config["riboseq_data"],
        project,
        sub_folder=sub_folder.as_posix(),
        note=note,
        is_unique=is_unique,
        fraction=fraction,
        reweighting_iterations=reweighting_iterations,
        is_filtered=is_filtered,
    )
    if args.overwrite or not Path(orfs_filen).is_file():
        bed_utils.write_bed(orfs[ext_field_names], orfs_filen)
    else:
        msg = f"Output file {orfs_filen} exists, skipping call!"
        logger.warning(msg)

    # Preparing Circos output
    msg = "Preparing Circos data"
    logger.info(msg)
    circos_filen = get_circos_graph(orfs, sub_folder, config, args)

    # preparing the BED track
    # add ORF colors
    msg = "Preparing IGV ORF track"
    logger.info(msg)

    # add color before, then cut to BED12
    for color, label in orf_type_colors.items():
        label_m = orfs["orf_category"] == label
        orfs.loc[label_m, "color"] = color
    orfs = orfs[bed_utils.bed12_field_names]
    orfs.drop_duplicates(inplace=True)

    igv_filen = filenames.get_riboseq_predicted_orfs(
        config["riboseq_data"],
        project,
        sub_folder=sub_folder.as_posix(),
        note=note,
        is_unique=is_unique,
        fraction=fraction,
        reweighting_iterations=reweighting_iterations,
        is_filtered=is_filtered,
    )
    igv_filen = igv_filen.replace("predicted-orfs", "igv-orfs")
    if args.overwrite or not Path(igv_filen).is_file():
        bed_utils.write_bed(orfs, igv_filen)
    else:
        msg = f"Output file {igv_filen} exists, skipping call!"
        logger.warning(msg)

    # log options for the dashboard
    dashboard_data = {
        "is_filtered": is_filtered,
        "min_samples": args.min_samples,
        "keep_other": args.keep_other,
        "no_replicates": args.no_replicates,
        "circos_bin_width": args.circos_bin_width,
        "circos_graph": circos_filen,
        "predicted_orfs": orfs_filen,
        "igv_tracks": {
            "orfs": {"file": igv_filen, "format": "bed"},
            "annotation": {"file": transcript_bed, "format": "bed"},
        },
        "date_time": "{:%Y-%m-%d %H:%M:%S}".format(datetime.datetime.now()),
    }
    filen = Path(sub_folder, f"{project}.summarize_options.json")
    filen = Path(config["riboseq_data"], filen)
    if args.overwrite or not filen.is_file():
        with open(filen, "w") as f:
            json.dump(dashboard_data, f)
    else:
        msg = f"Output file {filen} exists, skipping call!"
        logger.warning(msg)

    if args.show_orf_periodicity:
        msg = "Adding ORF periodicity figures"
        logger.info(msg)
        create_all_figures(config, sample_name_map, condition_name_map, args)


if __name__ == "__main__":
    main()
