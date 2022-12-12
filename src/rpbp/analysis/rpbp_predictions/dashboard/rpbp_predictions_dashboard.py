#! /usr/bin/env python3

import os
import argparse
import yaml
import json
import ast
import base64
import dash
import gzip

from pathlib import Path
from io import BytesIO
from dash import Dash, html, dcc, Input, Output, State, dash_table

import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.express as px

import pbiotools.ribo.ribo_filenames as filenames
import pbiotools.ribo.ribo_utils as ribo_utils

from rpbp.defaults import metagene_options, orf_type_colors, orf_type_labels, orf_type_name_map

#from controls import STACK_CTS_ORDER, STACK_CTS_NAME, FUNNEL_CTS_NAME

import dash_bio


import matplotlib.pyplot as plt
import seaborn as sns

sns.set({"ytick.direction": u'out'}, style = 'ticks') 
sns.set(rc = {'axes.facecolor': '#F9F9F8', 'figure.facecolor': '#F9F9F8'})


# ------------------------------------------------------ Functions ------------------------------------------------------


def parse_args():
    parser = argparse.ArgumentParser(description="""Launch a Dash app to
                                     visualize ORF predictions from Rp-Bp.""")

    parser.add_argument("config", type=str, 
                        help="Configuration (yaml) file used by Rp-Bp.")

    parser.add_argument("--debug", "-d", action="store_true", 
                        help="Enable debug mode")

    args = parser.parse_args()

    return args.config, args.debug

# TODO: do we need this?
# taken from https://github.com/4QuantOSS/DashIntro/blob/master/notebooks/Tutorial.ipynb
def fig_to_uri(in_fig, close_all=True, **save_args):
    # type: (plt.Figure) -> str
    """
    Save a figure as a URI
    :param in_fig:
    :return:
    """
    out_img = BytesIO()
    in_fig.savefig(out_img, format='png', **save_args)
    if close_all:
        in_fig.clf()
        plt.close('all')
    out_img.seek(0)  # rewind file
    encoded = base64.b64encode(out_img.read()).decode("ascii").replace("\n", "")
    return "data:image/png;base64,{}".format(encoded)


def fmt_tooltip(row):
    
    conditions = row.condition.split("|")
    bayes_factor_mean = row.bayes_factor_mean.split("|")
    bayes_factor_var = row.bayes_factor_var.split("|")
    profile_sum = row.profile_sum.split("|")
    in_frame = row.in_frame.split("|")
    
    it = zip(conditions,
       bayes_factor_mean,
       bayes_factor_var,
       profile_sum,
       in_frame
    )
    
    #fmt = []
    #for cond, bfm, bfv, ps, inf in it:
        #fmt.append(f"*{cond}*\nBF mean: {bfm}, BF var: {bfv}\n" \
                   #f"P-sites: {ps} ({inf}% in-frame)")
    #return "\n\n".join(fmt)
    
    fmt = ["| Condition | BF mean | BF var | P-sites | in-frame |"]
    fmt.append("| :-------- | :-----: | :----: | :-----: | -------: |")
    for cond, bfm, bfv, ps, inf in it:
        fmt.append(f"| {cond} | {bfm} | {bfv} | {ps} | {inf}% |")
    return "\n".join(fmt)


# ------------------------------------------------------ Set-up ------------------------------------------------------

# *** default path to summary data created by `summarize-rpbp-predictions`
sub_folder = Path("analysis", "rpbp_predictions")
igv_folder = Path("analysis", "igv")

# *** load configuration
configf, debug = parse_args()
config = yaml.load(open(configf), Loader=yaml.FullLoader)

project_name = config.get("project_name", "rpbp")
path_to_data = config["riboseq_data"]
path_to_genome = config["genome_base_path"]

prj_md_text = f"**Project name:** {project_name}\n\n" \
              f"**Path to genome:** {path_to_genome}\n\n" \
              f"**Path to data:** {path_to_data}\n\n" 

is_unique = not ('keep_riboseq_multimappers' in config)
config_note = config.get("note", None)

fraction = config.get('smoothing_fraction', None)
reweighting_iterations = config.get('smoothing_reweighting_iterations', None)
    
# pbiotools.ribo.ribo_utils._return_key_dict 
sample_name_map = ribo_utils.get_sample_name_map(config) # default to riboseq_samples.keys()
condition_name_map = ribo_utils.get_condition_name_map(config) # default to riboseq_biological_replicates.keys()
# work with pretty names, but we need the original names to retrieve files (e.g. metagene profiles)
reverse_sample_name_map = {sample_name_map[key]: key 
                           for key in config["riboseq_samples"].keys()}

# here some options, like metagene, etc.


# *** color palettes
pxcols = px.colors.qualitative.Set3
pal_frames = [pxcols[4], pxcols[3], pxcols[0]]
# TODO: pxcols to RGB 
pal_frames_meta = sns.palettes.color_palette(palette="Set3", n_colors=6)
pal_frames_meta = [pal_frames_meta[4], pal_frames_meta[3], pal_frames_meta[0]]
pal_bars = [pxcols[0], pxcols[1], pxcols[2]]

## *** app components


# ***
   


col_rev = {v:k for k,v in orf_type_colors.items()}
row_col = {}
for orf_type, labels in orf_type_labels.items():
    types = [orf_type_name_map[label] for label in labels]
    for t in types:
        row_col[t] = col_rev[orf_type]
        
style_data_conditional = [ {
    'if': {
        'filter_query': f'{{orf_type}} = "{t}"',
        'column_id': 'orf_type'
    },
    'backgroundColor': c,
    'color': 'white'
} for t, c in row_col.items()]
    
style_data_conditional.append({'if': {'column_id': 'transcripts'}, 'textOverflow': 'ellipsis', 'overflow': 'hidden', 'maxWidth': 0})
        
# *** load data

# 1. ORFs
filen = filenames.get_riboseq_predicted_orfs(
        config["riboseq_data"], 
        project_name, 
        sub_folder=sub_folder.as_posix(), 
        note=config_note,
        is_unique=is_unique,         
        fraction=fraction,
        reweighting_iterations=reweighting_iterations,
        is_filtered=True # TODO
    )
orfs = pd.read_csv(filen, sep="\t") # bed_utils
orfs.columns = orfs.columns.str.replace("#", "")
# TODO: prep unique table with information for tooltips, and list of transcripts in tooltip for long list

#orfs_columns = ['id', 'seqname', 'orf_len', 'orf_type'] 
#orfs = orfs[orfs_columns]

# AA - TODO check
orfs["orf_len"] = orfs["orf_len"]/3

orfs["profile_sum"] = orfs[["x_1_sum", "x_2_sum", "x_3_sum"]].sum(axis=1)
orfs["profile_sum"] = orfs["profile_sum"].astype(int)
orfs["in_frame"] = orfs["x_1_sum"].div(orfs["profile_sum"].values)*100
orfs["in_frame"] = orfs["in_frame"].apply(np.round).astype(int)

orfs["bayes_factor_mean"] = orfs["bayes_factor_mean"].apply(np.round).astype(int)
orfs["bayes_factor_var"] = orfs["bayes_factor_var"].apply(np.round).astype(int)

display_table = orfs.groupby("id", as_index=False)["condition"].agg({"condition": lambda x: "|".join(x)})

df = orfs.groupby("id", as_index=False)["bayes_factor_mean"].agg({"bayes_factor_mean": lambda x: "|".join([str(y) for y in x])})
display_table = display_table.join(df["bayes_factor_mean"])
df = orfs.groupby("id", as_index=False)["bayes_factor_var"].agg({"bayes_factor_var": lambda x: "|".join([str(y) for y in x])})
display_table = display_table.join(df["bayes_factor_var"])
df = orfs.groupby("id", as_index=False)["profile_sum"].agg({"profile_sum": lambda x: "|".join([str(y) for y in x])})
display_table = display_table.join(df["profile_sum"])
df = orfs.groupby("id", as_index=False)["in_frame"].agg({"in_frame": lambda x: "|".join([str(y) for y in x])})
display_table = display_table.join(df["in_frame"])


display_table["orf_info"] = display_table.apply(fmt_tooltip, axis=1)

DISPLAY_FIELDS = ["seqname", "id", "orf_len", "orf_type", "biotype", "transcripts", "gene_id", "gene_name", "gene_biotype"]


display_table = pd.merge(display_table, orfs[DISPLAY_FIELDS], on="id", how="left")
display_table = display_table[DISPLAY_FIELDS + ["orf_info"]]
#display_table = display_table[DISPLAY_FIELDS] # try now before adding tooltip
display_table.drop_duplicates(inplace=True) # after merge, we get duplicates...why?


sunburst_table = display_table.copy()
sunburst_table["length"] = "ORF"
sunburst_table.loc[sunburst_table["orf_len"]<100, "length"] = "sORF"
sunburst_table["count"] = 1
sunburst_col = row_col.copy()
sunburst_col["(?)"] = "#ededed"
sunburst_orfs = px.sunburst(
    sunburst_table,
    path=["length", "orf_type", "biotype"],
    values="count",
    color="orf_type",
    color_discrete_map=sunburst_col,
).update_traces(hovertemplate="%{label}<br>" + "Count: %{value}").update_layout(
    {
        "margin": dict(t=0, l=0, r=0, b=10),
        "plot_bgcolor": "rgba(0,0,0,0)", 
        "paper_bgcolor": "rgba(0,0,0,0)"
    }
)


#print(display_table)

filen = Path(path_to_data, igv_folder, f"{Path(config['fasta']).name}.txt")
with open(filen, 'r') as f:
    fastaURL = f.read()
    
filen = Path(path_to_data, igv_folder, f"{Path(config['fasta']).name}.fai.txt")
with open(filen, 'r') as f:
    indexURL = f.read()
    
filen = Path(path_to_data, igv_folder, f"{Path(config['gtf']).name}.txt")
with open(filen, 'r') as f:
    GTF = f.read()
    
filen = filenames.get_riboseq_predicted_orfs(
        config["riboseq_data"], 
        project_name, 
        sub_folder=igv_folder.as_posix(), 
        note=config_note,
        is_unique=is_unique,         
        fraction=fraction,
        reweighting_iterations=reweighting_iterations,
        is_filtered=True # TODO
    )
filen = filen.replace(".gz", ".txt")
with open(filen, 'r') as f:
    BED = f.read()    


# we need to know which ORF types are avavailable        
#option_orf_types = display_table.orf_type.unique()

option_orf_types = [
    {"label": x, "value": x} for x in display_table.orf_type.unique()
]

orf_type_default = "CDS" if any(["CDS" in d["label"] for d in option_orf_types]) else option_orf_types[0]['value']
drop_orf_types = dcc.Dropdown(
    id="drop_orf_types",
    clearable=False,
    searchable=False,
    options=option_orf_types,
    value=orf_type_default,
    style={"margin-top": "4px", "box-shadow": "0px 0px #73a5c8", "border-color": "#73a5c8"},
)

        

# we cannot do that on this one because it has multiple times same orf (per sample/conditions)
# orfs.sort_values(by=["orf_num"], inplace=True)
# orfs.set_index("orf_num", inplace=True, drop=True, verify_integrity=True)

# all in all, this will NOT be this df that we show in the table, but the unique one...
# if we have numpy we could use np.ceil
import math
PAGE_SIZE = 10
page_count = math.ceil(len(display_table)/PAGE_SIZE)



# IGV

_COMPONENT_ID = 'igv-chart'
#_COMPONENT_ID = 'default-igv'
    
reference={
        'id': 'IGV reference',
        'name': f'config["genome_name"]',
        'fastaURL': fastaURL,
        'indexURL': indexURL,
        'order': 1000000,
        'tracks': [
            {
                'name': 'Annotations',
                'url': GTF,
                'displayMode': 'EXPANDED',
                'nameField': 'gene',
                'height': 150,
                'color': 'rgb(0,0,0)',
                "type": "annotation",
                "format": "gff"
            },
            {
                'name': 'ORFs',
                'url': BED,
                'displayMode': 'EXPANDED',
                #'nameField': 'gene',
                'height': 150,
                "type": "annotation",
                "format": "bed"
            },
            
        ]
    }


filen = Path(path_to_data, igv_folder, f"{config['genome_name']}.circos_graph_data.json")
circos_graph_data = json.load(open(filen, "r"))


circos_layout_config = {
    "innerRadius": 150,
    "outerRadius": 200,
    "cornerRadius": 4,
    "labels": {
        "size": 15,
        "radialOffset": 200,
        "color": "#000000",
    },
    #"ticks": {
        #"color": "#000000",
        #"labelColor": "#000000",
        #"spacing": 10000000,
        #"labelSuffix": "Mb",
        #"labelDenominator": 1000000,
        #"labelSize": 10,
    #},
    "ticks": {"display": False},
}
circos_innerRadius = 1
circos_outerRadius = 2
circos_tracks_config = {"innerRadius": circos_innerRadius, 
                        "outerRadius": circos_outerRadius, 
                        "color": row_col[orf_type_default]}

# ------------------------------------------------------ APP ------------------------------------------------------

app = dash.Dash(__name__)

# do we have to expose the flask variable in the file?
# server = app.server

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
                # TODO: short intro, ref to Rp-Bp docs, etc.
                html.H1(
                    "Ribo-seq quality control",
                    style={"color": "rgb(0 0 0)"},
                ),
                html.Br(),
                dcc.Markdown(prj_md_text),
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
                                            id='igv-tabs',
                                            value='summary',
                                            children=[
                                                dcc.Tab(
                                                    label='Data summary',
                                                    value='summary',
                                                    children=html.Div(className='control-tab', children=[
                                                        dcc.Markdown(
                                                            """
                                                            Summary here
                                                            """
                                                        )
                                                    ])
                                                ),
                                                dcc.Tab(
                                                    label='ORF labels',
                                                    value='schematic',
                                                    children=html.Div(className='control-tab', children=[
                                                        dcc.Markdown(
                                                            """
                                                            OLF labels
                                                            """
                                                        )
                                                    ]),
                                                ),
                                            ],
                                        ),
                                    ],
                                    className="box"
                                ),  
                            ],
                            ),
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.H2("""1. ORF predictions per length, type, and host transcript biotype"""
                                        ),
                                        html.Br(), # ad hoc
                                        html.Br(),
                                        html.Br(),
                                        dcc.Graph(figure=sunburst_orfs),
                                    ],
                                    className="box",
                                    style={"width": "50%"},
                                ),
                                html.Div(
                                    [
                                        html.H2("""2. ORF type distribution along genomic coordinates"""
                                        ),
                                        html.Div(
                                            [
                                                html.Label("Select ORF type",
                                                        style={"margin": "10px"},
                                                ),
                                                html.Div(
                                                    [
                                                        drop_orf_types,
                                                    ],
                                                    style={
                                                        "width": "25%",
                                                        "margin-right": "25px"
                                                    },
                                                ),
                                            ],
                                            className="row",
                                        ),
                                        dash_bio.Circos(
                                            id="circos_fig",
                                            layout=circos_graph_data["genome"],
                                            config=circos_layout_config,
                                            tracks=[
                                                {
                                                    "type": "HISTOGRAM",
                                                    "data": circos_graph_data[f"histogram_{orf_type_default}"],
                                                    "config": circos_tracks_config,
                                                }
                                            ],
                                        ),
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
                                        html.H2("""3. ORF predictions (Table)"""
                                        ),
                                        # all backend paging/sorting/filtering
                                        dash_table.DataTable(
                                            id='datatable',
                                            columns=[
                                                {"name": i, "id": i} for i in DISPLAY_FIELDS
                                            ],
                                            page_current=0,
                                            page_size=PAGE_SIZE,
                                            page_count=page_count,
                                            page_action='custom',
                                            # sorting
                                            sort_action='custom',
                                            sort_mode='multi', # allow multi-column sorting
                                            sort_by=[],
                                            tooltip_header={'orf_len': 'Length in AA'},
                                            # Style headers with a dotted underline to indicate a tooltip
                                            # conditionally?
                                            #style_header={
                                                #'textDecoration': 'underline',
                                                #'textDecorationStyle': 'dotted',
                                            #},

                                            #tooltip_data=[{
                                                #'id': {'value': row['orf_info'], 'type': 'markdown'},
                                                #'transcripts': {'value': row['transcripts'], 'type': 'markdown'},
                                                #} for row in display_table.to_dict('records')
                                            #],
                                            #tooltip_data=tooltip,
                                            tooltip_delay=0,
                                            tooltip_duration=None,
                                            # Overflow into ellipsis
                                            #style_cell={
                                                #'overflow': 'hidden',
                                                #'textOverflow': 'ellipsis',
                                                #'maxWidth': 0,
                                            #},
                                            style_data_conditional=style_data_conditional,
                                        ),
                                        html.Label("""Hint: Hover over ORF ids to see in which sample/replicates they were found, 
                                            with evidence from Bayes factor and P-site counts. All transcripts compatible with
                                            an ORF are listed.""",
                                            style={"font-style": "italic"}
                                        ),
                                    ],
                                    # TODO: CSS/class style for table
                                    # className="dbc"
                                ),
                            ],
                            className="box",
                        ),
                        html.Div(
                            [
                                html.H2("""4. ORF predictions (IGV genome browser)"""
                                        ),
                                dash_bio.Igv(
                                    id=_COMPONENT_ID,
                                    reference=reference
                                ),
                            ],
                            className="box"
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
[Output('datatable', 'data'),
 Output('datatable', 'tooltip_data')],
[Input('datatable', "page_current"),
Input('datatable', "page_size"),
Input('datatable', 'sort_by')])
def update_table(page_current, page_size, sort_by):
    
    #return orfs.iloc[
        #page_current*page_size:(page_current+ 1)*page_size
    #].to_dict('records')
    
    if len(sort_by):
        #df = orfs.sort_values(
            #sort_by[0]['column_id'],
            #ascending=sort_by[0]['direction'] == 'asc',
            #inplace=False
        #)
        df = display_table.sort_values(
            [col['column_id'] for col in sort_by],
            ascending=[
                col['direction'] == 'asc'
                for col in sort_by
            ],
            inplace=False
        )
    else:
        # No sort is applied
        df = display_table

    data = df.iloc[
        page_current*page_size:(page_current+ 1)*page_size
    ].to_dict('records')
    
    tooltip_data=[{
                'id': {'value': row['orf_info'], 'type': 'markdown'},
                'transcripts': {'value': str(row['transcripts']), 'type': 'markdown'},
                } for row in data
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
        "color": row_col[value]
    }
    current[0].update(
        data=circos_graph_data[f"histogram_{value}"], 
        type="HISTOGRAM",
        config=tracks_config
    )
    return current
     

if __name__ == "__main__":
    app.run_server(debug=debug)
