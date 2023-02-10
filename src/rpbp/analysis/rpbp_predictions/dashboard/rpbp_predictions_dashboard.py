#! /usr/bin/env python3

import os
import argparse
import yaml
import json
import dash
import dash_bio
import flask
import threading
import webbrowser

from pathlib import Path
from dash import Dash, html, dcc, Input, Output, State, dash_table, ctx

import numpy as np
import pandas as pd
import plotly.express as px

from rpbp.defaults import orf_type_colors, orf_type_labels, orf_type_name_map

# ------------------------------------------------------ Functions ------------------------------------------------------


def get_parser():
    parser = argparse.ArgumentParser(
        description="Launch a Dash app to visualize "
        "Ribo-seq ORF predicted with Rp-Bp."
    )

    parser.add_argument(
        "config",
        type=str,
        help="A YAML configuration file." "The same used to run the pipeline.",
    )

    parser.add_argument("--debug", "-d", action="store_true", help="Enable debug mode.")

    parser.add_argument("--host", type=str, default="localhost", help="Host.")

    parser.add_argument("--port", type=int, default=8050, help="Port number.")

    args = parser.parse_args()

    return args.config, args.debug, args.host, args.port


def fmt_tooltip(row):

    conditions = row.condition.split("|")
    bayes_factor_mean = row.bayes_factor_mean.split("|")
    bayes_factor_var = row.bayes_factor_var.split("|")
    profile_sum = row.profile_sum.split("|")
    in_frame = row.in_frame.split("|")

    it = zip(conditions, bayes_factor_mean, bayes_factor_var, profile_sum, in_frame)

    fmt = ["| Condition | BF mean | BF var | P-sites | in-frame |"]
    fmt.append("| :-------- | :-----: | :----: | :-----: | -------: |")
    for cond, bfm, bfv, ps, inf in it:
        fmt.append(f"| {cond} | {bfm} | {bfv} | {ps} | {inf}% |")
    return "\n".join(fmt)


def get_orf_type_counts(condition):

    orf_type_counts = condition.groupby(["orf_type", "strand"]).size()
    orf_type_counts = orf_type_counts.reset_index(name="count")

    return orf_type_counts


def filter_sort_table(filter_query, sort_by):

    filtering_expressions = filter_query.split(" && ")
    df = display_table.copy()

    # first apply filtering
    for filter_part in filtering_expressions:
        col_name, operator, filter_value = split_filter_part(filter_part)

        if operator in ("eq", "lt", "le", "gt", "ge"):
            # these operators match pandas series operator method names
            df = df.loc[getattr(df[col_name], operator)(filter_value)]
        elif operator == "contains":
            df = df.loc[df[col_name].str.contains(filter_value)]

    if len(sort_by):
        df = df.sort_values(
            [col["column_id"] for col in sort_by],
            ascending=[col["direction"] == "asc" for col in sort_by],
            inplace=False,
        )

    return df


# ------------------------------------------------------ Set-up ------------------------------------------------------

# *** default path to summary data created by `summarize-rpbp-predictions`
sub_folder = Path("analysis", "rpbp_predictions")

# *** load configuration
configf, debug, host, port = get_parser()
config = yaml.load(open(configf), Loader=yaml.FullLoader)

project_name = config.get("project_name", "rpbp")
path_to_data = config["riboseq_data"]

prj_md_text = (
    f"**Project name:** {project_name}\n\n"
    f"**Data location:** {path_to_data}\n\n"
    f"---"
)

# *** load extended configuration, etc.
filen = Path(path_to_data, sub_folder, f"{project_name}.summarize_options.json")
config.update(json.load(open(filen, "r")))
config["orfs"] = config["igv_tracks"]["orfs"]["file"]
config["annotation"] = config["igv_tracks"]["annotation"]["file"]

is_filtered = config["is_filtered"]
is_filtered_str = "filtered" if is_filtered else "unfiltered"
no_repl_str = "excluding" if config["no_replicates"] else "including"
min_samples_str = (
    f"ORFs were *filtered* to keep only those "
    f"predicted in at least {config['min_samples']} samples."
)
default_str = (
    "ORFs were *not* filtered based on the number of predictions per sample (default)."
)
min_samples_str = min_samples_str if config["min_samples"] > 1 else default_str

results_md_text = (
    f"Ribo-seq ORFs from *{config['date_time']}*.\n\n"
    f"Using ORFs from *{is_filtered_str}* predictions, *{no_repl_str}* merged replicates.\n\n"
    f"{min_samples_str}\n\n"
    f"The bin width for counting ORF predictions along chromosomes is *{config['circos_bin_width']}b*."
)

labels_md_text = """
    **CDS**: Canonical (annotated) coding sequence\n
    **altCDS**: Alternative CDS (*e.g.* N/C-terminus extension/truncation, alternatively spliced variants, *etc.*)\n
    **intORF**: Translation event within a CDS (in- or out-of-frame)\n
    **uORF/uoORF**: Translation event in the 5' untranslated region (UTR) of or partially overlapping an annotated protein-coding gene\n
    **dORF/doORF**: Translation event in the 3' untranslated region (UTR) of or partially overlapping an annotated protein-coding gene\n
    **ncORF**: Translation event in an RNA annotated as non-coding (lncRNA, pseudogene, *etc.*)\n
    **Novel**: Translation event inter- or intragenic (only when Rp-Bp is run with a *de novo* assembly)
    """

col_rev = {v: k for k, v in orf_type_colors.items()}
row_col = {}
for orf_type, labels in orf_type_labels.items():
    types = [orf_type_name_map[label] for label in labels]
    for t in types:
        row_col[t] = col_rev[orf_type]

# *** load/wrangle data
orfs = pd.read_csv(config["predicted_orfs"], sep="\t", low_memory=False)  # bed_utils
orfs.columns = orfs.columns.str.replace("#", "")
orfs["orf_len"] = orfs["orf_len"] / 3
orfs["profile_sum"] = orfs[["x_1_sum", "x_2_sum", "x_3_sum"]].sum(axis=1)
orfs["profile_sum"] = orfs["profile_sum"].astype(int)
orfs["in_frame"] = orfs["x_1_sum"].div(orfs["profile_sum"].values) * 100
orfs["in_frame"] = orfs["in_frame"].apply(np.round).astype(int)
orfs["bayes_factor_mean"] = orfs["bayes_factor_mean"].apply(np.round).astype(int)
orfs["bayes_factor_var"] = orfs["bayes_factor_var"].apply(np.round).astype(int)

# main table - REDEFINE columns!
TABLE_FIELDS = [
    "seqname",
    "id",
    "orf_len",
    "orf_type",
    "biotype",
    "transcripts",
    "gene_id",
    "gene_name",
    "gene_biotype",
]
DISPLAY_FIELDS = [
    "Chrom",
    "ORF ID",
    "ORF length",
    "Category",
    "Transcript biotype",
    "Transcripts",
    "Gene ID",
    "Gene name",
    "Gene biotype",
]

display_table = orfs.groupby("id", as_index=False)["condition"].agg(
    {"condition": lambda x: "|".join(x)}
)
df = orfs.groupby("id", as_index=False)["bayes_factor_mean"].agg(
    {"bayes_factor_mean": lambda x: "|".join([str(y) for y in x])}
)
display_table = display_table.join(df["bayes_factor_mean"])
df = orfs.groupby("id", as_index=False)["bayes_factor_var"].agg(
    {"bayes_factor_var": lambda x: "|".join([str(y) for y in x])}
)
display_table = display_table.join(df["bayes_factor_var"])
df = orfs.groupby("id", as_index=False)["profile_sum"].agg(
    {"profile_sum": lambda x: "|".join([str(y) for y in x])}
)
display_table = display_table.join(df["profile_sum"])
df = orfs.groupby("id", as_index=False)["in_frame"].agg(
    {"in_frame": lambda x: "|".join([str(y) for y in x])}
)
display_table = display_table.join(df["in_frame"])
display_table["orf_info"] = display_table.apply(fmt_tooltip, axis=1)
# ad hoc - if we have the standardized ORFs
if "PHASE I ORFs" in orfs.columns:
    TABLE_FIELDS.extend(["PHASE I ORFs", "SS ORFs"])
    DISPLAY_FIELDS.extend(["Phase I", "Single-study"])
    orfs[["PHASE I ORFs", "SS ORFs"]] = orfs[["PHASE I ORFs", "SS ORFs"]].fillna(
        value="NA"
    )
# missing GTF fields - typically with de novo
orfs[TABLE_FIELDS[4:]] = orfs[TABLE_FIELDS[4:]].fillna(value="NA")
display_table = pd.merge(display_table, orfs[TABLE_FIELDS], on="id", how="left")
display_table = display_table[TABLE_FIELDS + ["orf_info"]]
display_table.drop_duplicates(inplace=True)

# summary table (tab) - show available samples and ORF types
orf_tab = orfs.groupby("condition").apply(get_orf_type_counts)
orf_tab.reset_index(inplace=True)
orf_tab.drop(columns="level_1", inplace=True)
orf_tab = orf_tab.pivot(
    index="condition", columns=["orf_type", "strand"], values="count"
)
orf_tab.fillna(0, inplace=True)
# reorder columns
all_multi = [
    (orf_type, strand)
    for orf_type in orf_type_name_map.values()
    for strand in ["+", "-"]
]
all_multi = [c for c in all_multi if c in orf_tab.columns]
orf_tab = orf_tab[all_multi]

# sunburst - no callback
sunburst_table = display_table.copy()
sunburst_table["length"] = "ORF"
sunburst_table.loc[sunburst_table["orf_len"] < 100, "length"] = "sORF"
sunburst_table["count"] = 1
sunburst_col = row_col.copy()
sunburst_col["(?)"] = "#ededed"
sunburst_orfs = (
    px.sunburst(
        sunburst_table,
        path=["length", "orf_type", "biotype"],  # transcript biotype
        values="count",
        color="orf_type",
        color_discrete_map=sunburst_col,
    )
    .update_traces(hovertemplate="%{label}<br>" + "Count: %{value}")
    .update_layout(
        {
            "margin": dict(t=0, l=0, r=0, b=10),
            "plot_bgcolor": "rgba(0,0,0,0)",
            "paper_bgcolor": "rgba(0,0,0,0)",
        }
    )
)

# main table - DISPLAY_FIELDS
display_table.rename(
    columns={k: v for k, v in zip(TABLE_FIELDS, DISPLAY_FIELDS)}, inplace=True
)

style_data_conditional = [
    {
        "if": {"filter_query": f'{{Category}} = "{t}"', "column_id": "Category"},
        "backgroundColor": c,
        "color": "white",
    }
    for t, c in row_col.items()
]

style_data_conditional.append(
    {
        "if": {"column_id": "Transcripts"},
        "textOverflow": "ellipsis",
        "overflow": "hidden",
        "maxWidth": 0,
    }
)

style_header_conditional = [
    {
        "if": {"column_id": f"{t}_{s}", "header_index": 0},
        "backgroundColor": c,
        "color": "white",
    }
    for t, c in row_col.items()
    for s in ["+", "-"]
]


tooltip_header = {
    "Chrom": "Chromosome",
    "ORF ID": "Rp-Bp ORF ID: Transcript:Chrom:Start-End:Strand",
    "ORF length": "ORF length (AA)",
    "Category": "Rio-seq ORF category (type or label)",
    "Transcript biotype": "Biotype of host transcript",
    "Transcripts": "Other compatible transcripts",
    "Gene ID": "Host gene ID",
    "Gene name": "Host gene name",
    "Gene biotype": "Host gene biotype",
}

PAGE_SIZE = 10
page_count = np.ceil(len(display_table) / PAGE_SIZE)

# table filtering
# from https://dash.plotly.com/datatable/callbacks
operators = [
    ["ge ", ">="],
    ["le ", "<="],
    ["lt ", "<"],
    ["gt ", ">"],
    ["eq ", "="],
    ["contains "],
]


def split_filter_part(filter_part):
    for operator_type in operators:
        for operator in operator_type:
            if operator in filter_part:
                name_part, value_part = filter_part.split(operator, 1)
                name = name_part[name_part.find("{") + 1 : name_part.rfind("}")]

                value_part = value_part.strip()
                v0 = value_part[0]
                if v0 == value_part[-1] and v0 in ("'", '"', "`"):
                    value = value_part[1:-1].replace("\\" + v0, v0)
                else:
                    try:
                        value = float(value_part)
                    except ValueError:
                        value = value_part

                # word operators need spaces after them in the filter string,
                # but we don't want these later
                return name, operator_type[0].strip(), value

    return [None] * 3


# data-dependent app components
option_orf_types = [{"label": x, "value": x} for x in display_table.Category.unique()]

orf_type_default = (
    "CDS"
    if any(["CDS" in d["label"] for d in option_orf_types])
    else option_orf_types[0]["value"]
)

drop_orf_types = dcc.Dropdown(
    id="drop_orf_types",
    clearable=False,
    searchable=False,
    options=option_orf_types,
    value=orf_type_default,
    style={
        "margin-top": "4px",
        "box-shadow": "0px 0px #73a5c8",
        "border-color": "#73a5c8",
    },
)

# IGV
_COMPONENT_ID = "igv-chart"
# _COMPONENT_ID = 'default-igv'

reference = {
    "id": "IGV reference",
    "name": config["genome_name"],
    "fastaURL": "data/fasta",
    "indexURL": "data/fasta/fai",
    "tracks": [
        {
            "name": "Ribo-seq ORFs",
            "format": config["igv_tracks"]["orfs"]["format"],
            "url": "data/orfs",
            "visibilityWindow": -1,
            "supportsWholeGenome": "false",
            "removable": "false",
            "order": 1000000,
        },
        {
            "name": "Annotation",
            "format": config["igv_tracks"]["annotation"]["format"],
            "url": "data/annotation",
            "visibilityWindow": -1,
            "supportsWholeGenome": "false",
            "removable": "false",
            "order": 1000000,
        },
    ],
}

# Circos
circos_graph_data = json.load(open(config["circos_graph"], "r"))

circos_layout_config = {
    "innerRadius": 150,
    "outerRadius": 200,
    "cornerRadius": 4,
    "labels": {
        "size": 15,
        "radialOffset": 200,
        "color": "#000000",
    },
    # "ticks": {
    # "color": "#000000",
    # "labelColor": "#000000",
    # "spacing": 10000000,
    # "labelSuffix": "Mb",
    # "labelDenominator": 1000000,
    # "labelSize": 10,
    # },
    "ticks": {"display": False},
}


circos_innerRadius = 1
circos_outerRadius = 2
circos_tracks_config = {
    "innerRadius": circos_innerRadius,
    "outerRadius": circos_outerRadius,
    "color": row_col[orf_type_default],
}  # "tooltipContent": {"name": "all"}

# ------------------------------------------------------ APP ------------------------------------------------------

app = dash.Dash(__name__)


@app.server.route("/data/<configval>", defaults={"suffix": None})
@app.server.route("/data/<configval>/<suffix>")
def config_data(configval, suffix):
    """Serve the file specified for the given key in the configuration.
    Potentially apply a suffix for derived files like an index.
    """

    # Extract the filename from configuration
    filename = config[configval]
    if suffix:
        filename = filename + "." + suffix

    return flask.send_file(filename)


app.layout = html.Div(
    [
        html.Div(
            [
                html.Img(
                    src=app.get_asset_url("logo-rpbp.png"),
                    style={
                        "position": "relative",
                        "width": "100%",
                        "left": "0px",
                        "top": "15px",
                    },
                ),
                html.H1(
                    "Ribo-seq ORF predictions",
                    style={"color": "rgb(0 0 0)"},
                ),
                html.Br(),
                dcc.Markdown(prj_md_text),
                dcc.Markdown(results_md_text),
            ],
            className="side_bar",
        ),
        html.Div(
            [
                html.Div(
                    [
                        html.Div(
                            [
                                html.Div(
                                    [
                                        dcc.Tabs(
                                            [
                                                dcc.Tab(
                                                    html.Div(
                                                        [
                                                            html.Div(
                                                                [
                                                                    html.H2(
                                                                        """ORF counts"""
                                                                    ),
                                                                    # table not interactive
                                                                    dash_table.DataTable(
                                                                        id="orf_tab",
                                                                        columns=[
                                                                            {
                                                                                "name": [
                                                                                    "Category",
                                                                                    "Strand",
                                                                                ],
                                                                                "id": "condition",
                                                                            }
                                                                        ]
                                                                        + [
                                                                            {
                                                                                "name": [
                                                                                    x1,
                                                                                    x2,
                                                                                ],
                                                                                "id": f"{x1}_{x2}",
                                                                            }
                                                                            for x1, x2 in orf_tab.columns
                                                                        ],
                                                                        data=[
                                                                            {
                                                                                **{
                                                                                    "condition": orf_tab.index[
                                                                                        n
                                                                                    ]
                                                                                },
                                                                                **{
                                                                                    f"{x1}_{x2}": y
                                                                                    for (
                                                                                        x1,
                                                                                        x2,
                                                                                    ), y in data
                                                                                },
                                                                            }
                                                                            for (
                                                                                n,
                                                                                data,
                                                                            ) in [
                                                                                *enumerate(
                                                                                    [
                                                                                        list(
                                                                                            x.items()
                                                                                        )
                                                                                        for x in orf_tab.T.to_dict().values()
                                                                                    ]
                                                                                )
                                                                            ]
                                                                        ],
                                                                        merge_duplicate_headers=True,
                                                                        style_header_conditional=style_header_conditional,
                                                                        style_data={
                                                                            "width": "100px",
                                                                            "maxWidth": "100px",
                                                                            "minWidth": "100px",
                                                                        },
                                                                        style_cell_conditional=[
                                                                            {
                                                                                "if": {
                                                                                    "column_id": "condition"
                                                                                },
                                                                                "width": "250px",
                                                                            },
                                                                        ],
                                                                        style_table={
                                                                            "overflowX": "auto"
                                                                        },
                                                                    ),
                                                                ],
                                                                style={
                                                                    "margin": "20px"
                                                                },
                                                            ),
                                                        ],
                                                        className="column",
                                                    ),
                                                    label="Data summary",
                                                ),
                                                dcc.Tab(
                                                    html.Div(
                                                        [
                                                            html.Div(
                                                                [
                                                                    html.Div(
                                                                        [
                                                                            html.Img(
                                                                                src=app.get_asset_url(
                                                                                    "schematic.png"
                                                                                ),
                                                                                style={
                                                                                    "width": "90%",
                                                                                    "position": "relative",
                                                                                    "margin": "15px",
                                                                                },
                                                                            ),
                                                                        ]
                                                                    ),
                                                                ],
                                                            ),
                                                            html.Div(
                                                                [
                                                                    dcc.Markdown(
                                                                        labels_md_text
                                                                    ),
                                                                ],
                                                                style={
                                                                    "width": "60%",
                                                                    "position": "relative",
                                                                    "margin": "5px",
                                                                    "margin-top": "15px",
                                                                },
                                                            ),
                                                        ],
                                                        className="row",
                                                    ),
                                                    label="Terminology and categories of Ribo-seq ORFs",
                                                ),
                                            ],
                                        ),
                                    ],
                                    className="box",
                                ),
                            ],
                        ),
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.H2(
                                            """1. ORF predictions per length, category, and host transcript biotype"""
                                        ),
                                        html.Br(),  # ad hoc
                                        html.Br(),
                                        dcc.Graph(
                                            figure=sunburst_orfs,
                                            style={
                                                "margin-top": "50px",
                                                "margin-bottom": "200px",
                                            },
                                        ),
                                        html.Label(
                                            """Hint: Click on sections to expand. Short ORFs (sORFs), also known
                                            as small ORFs (smORFs) are Ribo-seq ORFs < 100 amino acids in size. A transcript
                                            biotype is shown for the assigned host-transcript. Categories of Ribo-seq ORFs
                                            are assigned based on transcript-exon structure. To resolve seemingly incoherent
                                            assignments, look at all compatible transcripts in the ORF predictions table below.""",
                                            style={"font-style": "italic"},
                                        ),
                                    ],
                                    className="box",
                                    style={"width": "50%"},
                                ),
                                html.Div(
                                    [
                                        html.H2(
                                            """2. Distribution of Ribo-seq ORFs along genomic coordinates"""
                                        ),
                                        html.Div(
                                            [
                                                html.Label(
                                                    "Select category",
                                                    style={"margin": "10px"},
                                                ),
                                                html.Div(
                                                    [
                                                        drop_orf_types,
                                                    ],
                                                    style={
                                                        "width": "25%",
                                                        "margin-right": "25px",
                                                    },
                                                ),
                                            ],
                                            className="row",
                                        ),
                                        dash_bio.Circos(
                                            id="circos_fig",
                                            layout=circos_graph_data["genome"],
                                            config=circos_layout_config,
                                            # selectEvent={"0": "hover"},
                                            tracks=[
                                                {
                                                    "data": circos_graph_data[
                                                        f"histogram_{orf_type_default}"
                                                    ],
                                                    "type": "HISTOGRAM",
                                                    "config": circos_tracks_config,
                                                }
                                            ],
                                        ),
                                        # html.Div(id="circos_output"),
                                    ],
                                    className="box",
                                    style={"width": "50%"},
                                ),
                            ],
                            className="row",
                        ),
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.H2("""3. ORF predictions (Table)"""),
                                        html.Div(
                                            [
                                                html.Button(
                                                    "Download CSV", id="btn_csv"
                                                ),
                                                dcc.Download(
                                                    id="download-dataframe-csv"
                                                ),
                                            ]
                                        ),
                                        html.Br(),
                                        # all backend paging/sorting/filtering
                                        dash_table.DataTable(
                                            id="main_datatable",
                                            columns=[
                                                {"name": i, "id": i}
                                                for i in DISPLAY_FIELDS
                                            ],
                                            page_current=0,
                                            page_size=PAGE_SIZE,
                                            page_count=page_count,
                                            page_action="custom",
                                            # sorting
                                            sort_action="custom",
                                            sort_mode="multi",  # allow multi-column sorting
                                            sort_by=[],
                                            # filtering
                                            filter_action="custom",
                                            filter_query="",
                                            tooltip_header=tooltip_header,
                                            tooltip_delay=0,
                                            tooltip_duration=None,
                                            style_header={
                                                "textDecoration": "underline",
                                                "textDecorationStyle": "dotted",
                                            },
                                            style_data_conditional=style_data_conditional,
                                            # export_format="csv",
                                        ),
                                        dcc.Markdown(
                                            """**Filtering syntax:** Use *eq (=)*, *le (<=)*, *lt (<)*, *ge (>=)*
                                            , *gt (>)*, or *contains*, such as `< 100` in the "ORF length" column, or
                                            `contains C30F8.2.1` in the "Transcripts" column, or `= uORF` in the "Category"
                                            column, *etc*. The default filtering behavior depends on the data type *e.g.*
                                            `contains 12652389` in "ORF ID" will not work, but `contains 12652389-12654351`
                                            will. Check the little *Aa* box to toggle case sensitivity. Filters and sorting
                                            can be combined, but do not forget to clear the filters, this is not done
                                            automatically! You can download the full table, or a selection based on filters."""
                                        ),
                                        dcc.Markdown(
                                            """**Standardized Ribo-seq ORFs:** If present, additional columns are shown. Search
                                            for PHASE I ORFs using *e.g* "c1" for ORFs on chromosome 1, or just "orf". Search
                                            for Single-study ORFs using the same syntax, or "norep"."""
                                        ),
                                        html.Br(),
                                        html.Label(
                                            """Hint: Hover over ORF IDs to see in which sample/replicates they were found,
                                            with evidence from Bayes factors and P-site counts. All transcripts compatible with
                                            an ORF are listed (Transcripts).""",
                                            style={"font-style": "italic"},
                                        ),
                                    ],
                                ),
                            ],
                            className="box",
                        ),
                        html.Div(
                            [
                                html.H2("""4. ORF predictions (IGV genome browser)"""),
                                dash_bio.Igv(id=_COMPONENT_ID, reference=reference),
                            ],
                            className="box",
                        ),
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.P(
                                            [
                                                html.A(
                                                    "Dieterich Lab",
                                                    href="https://github.com/dieterich-lab/",
                                                    target="_blank",
                                                ),
                                                html.Br(),
                                                "Etienne Boileau",
                                            ],
                                            style={"font-size": "12px"},
                                        ),
                                    ],
                                    style={"width": "60%"},
                                ),
                                html.Div(
                                    [
                                        html.P(
                                            [
                                                "Rp-Bp",
                                                html.Br(),
                                                html.A(
                                                    "Documentation",
                                                    href="https://rp-bp.readthedocs.io/en/latest/index.html",
                                                    target="_blank",
                                                ),
                                                ", ",
                                                html.A(
                                                    "Source code",
                                                    href="https://github.com/dieterich-lab/rp-bp/",
                                                    target="_blank",
                                                ),
                                                ", ",
                                                html.A(
                                                    "Publication",
                                                    href="https://academic.oup.com/nar/article/45/6/2960/2953491",
                                                    target="_blank",
                                                ),
                                            ],
                                            style={"font-size": "12px"},
                                        )
                                    ],
                                    style={"width": "37%"},
                                ),
                            ],
                            className="footer",
                            style={"display": "flex"},
                        ),
                    ],
                    className="main",
                ),
            ]
        ),
    ]
)

# ------------------------------------------------------ Callbacks ------------------------------------------------------


@app.callback(
    [Output("main_datatable", "data"), Output("main_datatable", "tooltip_data")],
    [
        Input("main_datatable", "page_current"),
        Input("main_datatable", "page_size"),
        Input("main_datatable", "sort_by"),
        Input("main_datatable", "filter_query"),
    ],
)
def update_main_table(page_current, page_size, sort_by, filter_query):
    def refmt(trx):
        trxs = trx.split(",")
        return "\n".join(list(map(" ".join, zip(trxs[::2], trxs[1::2]))))

    df = filter_sort_table(filter_query, sort_by)

    data = df.iloc[page_current * page_size : (page_current + 1) * page_size].to_dict(
        "records"
    )

    tooltip_data = [
        {
            "ORF ID": {"value": row["orf_info"], "type": "markdown"},
            "Transcripts": {"value": refmt(row["Transcripts"]), "type": "markdown"},
        }
        for row in data
    ]

    return data, tooltip_data


@app.callback(
    Output("circos_fig", "tracks"),
    Input("drop_orf_types", "value"),
    State("circos_fig", "tracks"),
)
def hist_orf_type(value, current):

    tracks_config = {
        "innerRadius": circos_innerRadius,
        "outerRadius": circos_outerRadius,
        "color": row_col[value],
    }
    current[0].update(
        data=circos_graph_data[f"histogram_{value}"],
        type="HISTOGRAM",
        config=tracks_config,
    )
    return current


# @app.callback(
# Output("circos_output", "children"),
# Input("circos_fig", "eventDatum"),
# )
# def update_output(value):
# if value is not None:
# return [html.Div("{}: {}".format(v.title(), value[v])) for v in value.keys()]
# return "Hover over a bar to get more information."


@app.callback(
    Output("download-dataframe-csv", "data"),
    [
        Input("btn_csv", "n_clicks"),
        Input("main_datatable", "sort_by"),
        Input("main_datatable", "filter_query"),
    ],
    # State("main_datatable", "data"),
    prevent_initial_call=True,
)
def func(n_clicks, sort_by, filter_query):  # table_data

    # df = pd.DataFrame.from_dict(table_data)
    changed_inputs = [x["prop_id"] for x in ctx.triggered]
    if "btn_csv.n_clicks" in changed_inputs:
        df = filter_sort_table(filter_query, sort_by)
        df.drop(columns="orf_info", inplace=True)
        return dcc.send_data_frame(df.to_csv, "selected-orfs.csv", index=False)


def main():
    if "DISPLAY" in os.environ and os.environ["DISPLAY"]:
        threading.Timer(
            1, lambda: webbrowser.open_new(f"http://{host}:{port}/")
        ).start()
    app.run(debug=debug, host=host, port=port)


if __name__ == "__main__":
    main()
