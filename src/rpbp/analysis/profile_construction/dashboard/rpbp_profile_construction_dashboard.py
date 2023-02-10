#! /usr/bin/env python3

import os
import sys

import argparse
import yaml
import json
import ast
import base64
import dash
import threading
import webbrowser

from pathlib import Path
from io import BytesIO
from dash import Dash, html, dcc, Input, Output

import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.express as px

import rpbp.ribo_utils.utils as ribo_utils
import rpbp.ribo_utils.filenames as filenames

from rpbp.defaults import metagene_options

from rpbp.analysis.profile_construction.dashboard.controls import (
    STACK_CTS_ORDER,
    STACK_CTS_NAME,
    FUNNEL_CTS_NAME,
)

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

sns.set({"ytick.direction": "out"}, style="ticks")
sns.set(rc={"axes.facecolor": "#F9F9F8", "figure.facecolor": "#F9F9F8"})


# ------------------------------------------------------ Functions ------------------------------------------------------


def get_parser():
    parser = argparse.ArgumentParser(
        description="Launch a Dash app for quality "
        "control and visualization of ribosome profiling "
        "data processed with Rp-Bp."
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


def get_diff_counts(data_np):

    # add an extra column so the diff counts will work
    zeros = np.zeros((data_np.shape[0], 1))
    data_np = np.append(zeros, data_np, axis=1)
    # get the diffs so the stacks work correctly
    diff = np.diff(data_np)

    return diff


def get_window_counts(
    metagene_profile,
    offset,
    start_upstream_window,
    start_downstream_window,
    end_upstream_window,
    end_downstream_window,
):

    # profile around start codon
    start_upstream = start_upstream_window - offset
    start_downstream = start_downstream_window - offset

    mask_starts = metagene_profile["type"] == "start"
    m_start_upstream = metagene_profile["position"] >= start_upstream
    m_start_downstream = metagene_profile["position"] <= start_downstream

    mask_starts = m_start_upstream & m_start_downstream & mask_starts
    start_counts = metagene_profile.loc[mask_starts, "count"].values

    # profile around stop codon
    end_upstream = end_upstream_window - offset
    end_downstream = end_downstream_window - offset

    mask_ends = metagene_profile["type"] == "end"
    m_end_upstream = metagene_profile["position"] >= end_upstream
    m_end_downstream = metagene_profile["position"] <= end_downstream

    mask_ends = m_end_upstream & m_end_downstream & mask_ends
    end_counts = metagene_profile.loc[mask_ends, "count"].values

    return start_counts, end_counts


def get_profiles_bars(
    samples,
    start_upstream_window,
    start_downstream_window,
    end_upstream_window,
    end_downstream_window,
):

    start_window_size = start_downstream_window - start_upstream_window + 1
    end_window_size = end_downstream_window - end_upstream_window + 1

    metagene_profile_start = np.zeros(start_window_size)
    metagene_profile_end = np.zeros(end_window_size)

    for sample in samples:

        sample_name = reverse_sample_name_map[sample]

        # TODO: if we have the df with all lengths and offsets ready, we don't need to
        # read them each time, but just subset the pre-loaded dataframe
        lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(
            config,
            sample_name,
            is_unique=is_unique,  # defaults same as rp-bp set in ribo_utils TODO default_params=metagene_options
        )

        metagene_profiles = pd.read_csv(
            filenames.get_metagene_profiles(
                config["riboseq_data"],
                sample_name,
                is_unique=is_unique,
                note=config.get("note", None),
            )
        )

        for length, offset in zip(lengths, offsets):

            m_length = metagene_profiles["length"] == int(length)
            metagene_profile = metagene_profiles[m_length]

            start_counts, end_counts = get_window_counts(
                metagene_profile,
                int(offset),
                start_upstream_window,
                start_downstream_window,
                end_upstream_window,
                end_downstream_window,
            )
            metagene_profile_start += start_counts
            metagene_profile_end += end_counts

    return metagene_profile_start, metagene_profile_end


def get_profiles_lines(
    samples,
    start_upstream_window,
    start_downstream_window,
    end_upstream_window,
    end_downstream_window,
):

    start_window_size = start_downstream_window - start_upstream_window + 1
    end_window_size = end_downstream_window - end_upstream_window + 1

    all_ribo_metagene_profile_start_df = pd.DataFrame(columns=range(start_window_size))
    all_ribo_metagene_profile_end_df = pd.DataFrame(columns=range(end_window_size))

    for sample in samples:

        sample_name = reverse_sample_name_map[sample]

        # TODO: if we have the df with all lengths and offsets ready, we don't need to
        # read them each time, but just subset the pre-loaded dataframe
        lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(
            config,
            sample_name,
            is_unique=is_unique,  # defaults same as rp-bp set in ribo_utils
        )

        metagene_profiles = pd.read_csv(
            filenames.get_metagene_profiles(
                config["riboseq_data"],
                sample_name,
                is_unique=is_unique,
                note=config.get("note", None),
            )
        )

        window_count_start, window_count_end = 0, 0
        metagene_profile_start = np.zeros(start_window_size)
        metagene_profile_end = np.zeros(end_window_size)

        for length, offset in zip(lengths, offsets):

            m_length = metagene_profiles["length"] == int(length)
            metagene_profile = metagene_profiles[m_length]

            start_counts, end_counts = get_window_counts(
                metagene_profile,
                int(offset),
                start_upstream_window,
                start_downstream_window,
                end_upstream_window,
                end_downstream_window,
            )

            metagene_profile_start += start_counts
            window_count_start += np.sum(start_counts)
            metagene_profile_end += end_counts
            window_count_end += np.sum(end_counts)

        # collect all profiles for this sample
        metagene_profile_start /= window_count_start
        metagene_profile_end /= window_count_end
        all_ribo_metagene_profile_start_df.loc[sample_name] = metagene_profile_start
        all_ribo_metagene_profile_end_df.loc[sample_name] = metagene_profile_end

    ribo_start = all_ribo_metagene_profile_start_df.median(axis=0)
    ribo_end = all_ribo_metagene_profile_end_df.median(axis=0)

    return ribo_start, ribo_end


def metagene_plot_bars(
    metagene_profile_start,
    metagene_profile_end,
    start_upstream_window,
    start_downstream_window,
    end_upstream_window,
    end_downstream_window,
    step,
    figsize,
    y_label,
    y_max=None,
    layout=None,
):

    indices = [(0, 1, 2), (1, 0, 2), (2, 1, 0)]
    start_window_size = start_downstream_window - start_upstream_window + 1
    end_window_size = end_downstream_window - end_upstream_window + 1

    if y_max is None:
        ymax = max(max(metagene_profile_start), max(metagene_profile_end))
    else:
        ymax = y_max
    ymax = 1.1 * ymax

    xticklabels_start = range(start_upstream_window, start_downstream_window + 1, step)
    xticks_start = range(0, start_window_size, step)

    xticklabels_end = range(end_upstream_window, end_downstream_window + 1, step)
    xticks_end = range(0, end_window_size, step)

    idx_start = indices[abs(start_upstream_window) % 3]
    idx_end = indices[abs(end_upstream_window) % 3]
    # fixed size
    fig, axes = plt.subplots(ncols=2, figsize=figsize, sharey=True, layout=layout)

    # first, the start counts
    x_1 = metagene_profile_start[idx_start[0] :: 3]
    x_2 = metagene_profile_start[idx_start[1] :: 3]
    x_3 = metagene_profile_start[idx_start[2] :: 3]

    x_1_pos = np.arange(len(metagene_profile_start))[idx_start[0] :: 3]
    x_2_pos = np.arange(len(metagene_profile_start))[idx_start[1] :: 3]
    x_3_pos = np.arange(len(metagene_profile_start))[idx_start[2] :: 3]

    x_1_rects = axes[0].bar(x_1_pos, x_1, width=1, color=pal_frames[0])
    x_2_rects = axes[0].bar(x_2_pos, x_2, width=1, color=pal_frames[1])
    x_3_rects = axes[0].bar(x_3_pos, x_3, width=1, color=pal_frames[2])

    # now, the end counts
    x_1 = metagene_profile_end[idx_end[0] :: 3]
    x_2 = metagene_profile_end[idx_end[1] :: 3]
    x_3 = metagene_profile_end[idx_end[2] :: 3]

    x_1_pos = np.arange(len(metagene_profile_end))[idx_end[0] :: 3]
    x_2_pos = np.arange(len(metagene_profile_end))[idx_end[1] :: 3]
    x_3_pos = np.arange(len(metagene_profile_end))[idx_end[2] :: 3]

    x_1_rects = axes[1].bar(x_1_pos, x_1, width=1, color=pal_frames[0])
    x_2_rects = axes[1].bar(x_2_pos, x_2, width=1, color=pal_frames[1])
    x_3_rects = axes[1].bar(x_3_pos, x_3, width=1, color=pal_frames[2])

    axes[0].set_xticks(xticks_start)
    axes[0].set_xticklabels(xticklabels_start)
    axes[0].yaxis.tick_left()
    axes[0].set_ylim((0, ymax))

    axes[0].set_ylabel(y_label)  # scaling?

    axes[1].set_xticks(xticks_end)
    axes[1].set_xticklabels(xticklabels_end)
    axes[1].axes.get_yaxis().set_visible(False)

    sns.despine(ax=axes[0])
    sns.despine(ax=axes[1], left=True)

    return (fig, axes)


def metagene_plot_lines(
    metagene_profile_start,
    metagene_profile_end,
    start_upstream_window,
    start_downstream_window,
    end_upstream_window,
    end_downstream_window,
    step,
    figsize,
    y_label,
    y_max=None,
):

    start_window_size = start_downstream_window - start_upstream_window + 1
    end_window_size = end_downstream_window - end_upstream_window + 1

    if y_max is None:
        ymax = max(max(metagene_profile_start), max(metagene_profile_end))
    else:
        ymax = y_max
    ymax = 1.1 * ymax

    xticklabels_start = range(start_upstream_window, start_downstream_window + 1, step)
    xticks_start = range(0, start_window_size, step)

    xticklabels_end = range(end_upstream_window, end_downstream_window + 1, step)
    xticks_end = range(0, end_window_size, step)

    # fixed size
    fig, axes = plt.subplots(ncols=2, figsize=figsize, sharey=True)

    # first, the start counts
    x_pos = np.arange(start_window_size)
    x_rects = axes[0].plot(x_pos, metagene_profile_start, color=pal_frames[0])

    # now, the end counts
    # stop site, ribo, then rna
    x_pos = np.arange(end_window_size)
    x_rects = axes[1].plot(x_pos, metagene_profile_end, color=pal_frames[0])

    axes[0].set_xticks(xticks_start)
    axes[0].set_xticklabels(xticklabels_start)
    axes[0].yaxis.tick_left()
    axes[0].set_ylim((0, ymax))

    axes[0].set_ylabel(y_label)  # scaling?

    axes[1].set_xticks(xticks_end)
    axes[1].set_xticklabels(xticklabels_end)
    axes[1].axes.get_yaxis().set_visible(False)

    return (fig, axes)


# taken from https://github.com/4QuantOSS/DashIntro/blob/master/notebooks/Tutorial.ipynb
def fig_to_uri(in_fig, close_all=True, **save_args):
    # type: (plt.Figure) -> str
    """
    Save a figure as a URI
    :param in_fig:
    :return:
    """
    out_img = BytesIO()
    in_fig.savefig(out_img, format="png", **save_args)
    if close_all:
        in_fig.clf()
        plt.close("all")
    out_img.seek(0)  # rewind file
    encoded = base64.b64encode(out_img.read()).decode("ascii").replace("\n", "")
    return "data:image/png;base64,{}".format(encoded)


# ------------------------------------------------------ Set-up ------------------------------------------------------

# *** default path to summary data created by `summarize-rpbp-profile-construction`
sub_folder = Path("analysis", "profile_construction")

# *** load configuration
configf, debug, host, port = get_parser()
config = yaml.load(open(configf), Loader=yaml.FullLoader)

project_name = config.get("project_name", "rpbp")
path_to_data = config["riboseq_data"]
path_to_genome = config["genome_base_path"]

prj_md_text = (
    f"**Project name:** {project_name}\n\n"
    f"**Genome location:** {path_to_genome}\n\n"
    f"**Data location:** {path_to_data}\n\n"
    f"---"
)

is_unique = not config.get("keep_riboseq_multimappers", False)
config_note = config.get("note", None)

# ribo_utils._return_key_dict
sample_name_map = ribo_utils.get_sample_name_map(
    config
)  # default to riboseq_samples.keys()
condition_name_map = ribo_utils.get_condition_name_map(
    config
)  # default to riboseq_biological_replicates.keys()
# work with pretty names, but we need the original names to retrieve files (e.g. metagene profiles)
reverse_sample_name_map = {
    sample_name_map[key]: key for key in config["riboseq_samples"].keys()
}

# TODO: can we get min/max from data interactively?
start_upstream = config.get(
    "metagene_start_upstream", metagene_options["metagene_start_upstream"]
)
start_downstream = config.get(
    "metagene_start_downstream",
    metagene_options["metagene_start_downstream"],
)
end_upstream = config.get(
    "metagene_end_upstream",
    metagene_options["metagene_end_upstream"],
)
end_downstream = config.get(
    "metagene_end_downstream",
    metagene_options["metagene_end_downstream"],
)
# options allow both to be different, but for plotting
# we ensure a symmetric window, so restrict choice to
# the minimum length
min_upstream = min(start_upstream, end_upstream)
min_downstream = min(start_downstream, end_downstream)

# for read length profiles, we fix the window size
# TODO: this may not work if e.g. inconsistent values
# were set in the config...
read_length_start_upstream = -50
read_length_start_downstream = 21
read_length_end_upstream = -50
read_length_end_downstream = 21
read_length_step = 10

# *** color palettes
pal_blind = sns.color_palette("colorblind").as_hex()
pal_blind = [
    pal_blind[0],
    pal_blind[3],
    pal_blind[2],
    pal_blind[8],
    pal_blind[4],
    pal_blind[1],
]
pal_frames = pal_blind[:3]
pal_bars = pal_blind[2:5]
pal_set3 = pal_bars + [pal_blind[1], pal_blind[0], pal_blind[5]]

# *** app components
option_unique = [dict(label="All mapped reads", value="False")]
if is_unique:
    option_unique.append(dict(label="Uniquely mapped reads", value="True"))

option_samples = [
    dict(
        label="All",
        value=",".join(
            [sample_name_map[sample] for sample in config["riboseq_samples"].keys()]
        ),
    )
]
if "riboseq_biological_replicates" in config:
    if config["riboseq_biological_replicates"] is not None:
        riboseq_replicates = config["riboseq_biological_replicates"].items()
        for condition, replicates in riboseq_replicates:
            option_samples.append(
                dict(
                    label=condition_name_map[condition],
                    value=",".join([sample_name_map[sample] for sample in replicates]),
                )
            )

drop_samples = dcc.Dropdown(
    id="drop_samples",
    clearable=False,
    searchable=False,
    options=option_samples,
    value=option_samples[0]["value"],
    style={
        "margin-top": "4px",
        "box-shadow": "0px 0px #73a5c8",
        "border-color": "#73a5c8",
    },
)

drop_unique = dcc.Dropdown(
    id="drop_unique",
    clearable=False,
    searchable=False,
    options=option_unique,
    value=option_unique[0]["value"],
    style={"margin": "4px", "box-shadow": "0px 0px #73a5c8", "border-color": "#73a5c8"},
)

radio_stacked_reads = dcc.RadioItems(
    id="radio_stacked_reads",
    # className="radio", # customize with css
    options=[
        dict(
            label="All reads", title="Show the number of reads after each step", value=0
        ),
        dict(
            label="Exclude ribosomal/poor quality",
            title="Show the number of reads after each step, excluding reads "
            "of poor quality and those that map to ribosomal sequences",
            value=1,
        ),
    ],
    value=0,
    inline=True,
)

dropdown_available_samples = dcc.Dropdown(
    id="dropdown_available_samples",
    clearable=False,
    searchable=True,
    style={"margin": "4px", "box-shadow": "0px 0px #73a5c8", "border-color": "#73a5c8"},
)

read_length_slider = dcc.RangeSlider(
    id="read_length_slider",
    step=1,
    marks=None,
    tooltip={"placement": "bottom", "always_visible": True},
)

toggle_log = dbc.Switch(
    id="toggle_log",
    # label="Log scale",
    value=False,
)

# *** load data
# 1. read filtering counts
filename = filenames.get_riboseq_read_filtering_counts(
    config["riboseq_data"],
    project_name,
    sub_folder=sub_folder.as_posix(),
    note=config_note,
)
read_filtering_counts = pd.read_csv(filename)
read_filtering_counts["Sample"] = read_filtering_counts["sample"].apply(
    lambda x: sample_name_map[x]
)
read_filtering_counts = read_filtering_counts.sort_values("Sample")

# 2. read length distributions
filename = filenames.get_riboseq_read_length_distribution(
    config["riboseq_data"],
    project_name,
    sub_folder=sub_folder.as_posix(),
    note=config_note,
)
read_length_distributions = pd.read_csv(filename)
read_length_distributions["sample"] = read_length_distributions["sample"].apply(
    lambda x: sample_name_map[x]
)
read_length_distributions = read_length_distributions.melt(
    id_vars=["sample", "is_unique"], var_name="length", value_name="count"
)
read_length_distributions.rename(columns={"sample": "Sample"}, inplace=True)
read_length_distributions["length"] = read_length_distributions["length"].astype(int)

# 3. periodic offsets
filename = filenames.get_periodic_offsets(
    config["riboseq_data"],
    project_name,
    sub_folder=sub_folder.as_posix(),
    is_unique=is_unique,
    note=config_note,
)
lengths_and_offsets = pd.read_csv(filename)
lengths_and_offsets["sample"] = lengths_and_offsets["sample"].apply(
    lambda x: sample_name_map[x]
)
lengths_and_offsets.rename(columns={"sample": "Sample"}, inplace=True)

# 4. frame counts
filename = filenames.get_riboseq_frame_counts(
    config["riboseq_data"],
    project_name,
    sub_folder=sub_folder.as_posix(),
    is_unique=is_unique,
    note=config_note,
)
frame_counts = pd.read_csv(filename)
frame_counts["sample"] = frame_counts["sample"].apply(lambda x: sample_name_map[x])
frame_counts = frame_counts.sort_values("sample")
frame_counts = pd.melt(
    frame_counts,
    id_vars="sample",
    value_vars=["frame", "frame+1", "frame+2"],
    var_name="Frame",
)
frame_counts["p_site"] = (
    frame_counts.groupby(["sample"])["value"].transform(lambda x: x / x.sum()) * 100
)
frame_map = {"frame": "In-frame", "frame+1": "In-frame + 1", "frame+2": "In-frame + 2"}
frame_counts.replace({"Frame": frame_map}, inplace=True)
frame_counts.rename(columns={"sample": "Sample"}, inplace=True)


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
                html.Br(),
                html.Label(
                    "Select samples (conditions or replicates) to use:",
                    style={"font-weight": "bold"},  # how to padding-bottom
                ),
                drop_samples,
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
                                        html.H2(
                                            """1. Stacked bars and funnel charts showing
                                            the number of reads after each step of the
                                            filtering pipeline"""
                                        ),
                                        html.Br(),
                                        radio_stacked_reads,
                                        dcc.Graph(
                                            id="stacked_reads_fig",
                                        ),
                                        html.Label(
                                            """Hint: Hover over bars to create a funnel
                                            chart (right) for a given sample.""",
                                            style={"font-style": "italic"},
                                        ),
                                    ],
                                    className="box",
                                    style={"width": "60%"},
                                ),
                                html.Div(
                                    [
                                        dcc.Graph(id="funnel_fig"),
                                        html.Br(),  # why does the label sticks to the plot?
                                        html.Br(),  # TODO: flush to bottom
                                        html.Br(),
                                        html.Br(),
                                        html.Br(),
                                        html.Br(),
                                        html.Br(),
                                        html.Label(
                                            """Hint: Hover over chart to see proportions.""",
                                            style={
                                                "font-style": "italic",
                                                "vertical-align": "bottom",
                                            },
                                        ),
                                    ],
                                    className="box",
                                    style={"width": "40%"},
                                ),
                            ],
                            className="row",
                        ),
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.H2(
                                            """2. Per sample read length distribution and
                                        metagene profiles"""
                                        ),
                                        html.Div(
                                            [
                                                html.Label(
                                                    "Select one sample",
                                                    style={"margin": "10px"},
                                                ),
                                                html.Div(
                                                    [
                                                        dropdown_available_samples,
                                                    ],
                                                    style={
                                                        "width": "25%",
                                                        "margin-right": "25px",
                                                    },
                                                ),
                                                html.Label(
                                                    "Select alignments",
                                                    style={"margin": "10px"},
                                                ),
                                                html.Div(
                                                    [
                                                        drop_unique,
                                                    ],
                                                    style={"width": "25%"},
                                                ),
                                            ],
                                            className="row",
                                        ),
                                        html.Br(),
                                        html.Div(
                                            [
                                                html.Label(
                                                    "Log scale",
                                                    style={"margin-left": "10px"},
                                                ),
                                                html.Div(
                                                    [toggle_log],
                                                    style={"width": "10%"},
                                                ),
                                                html.Label("Select range"),
                                                html.Div(
                                                    [read_length_slider],
                                                    style={"width": "80%"},
                                                ),
                                            ],
                                            className="row",
                                        ),
                                        dcc.Graph(id="length_bars_fig"),
                                        html.Label(
                                            """Hint: Hover over bars to see the metagene
                                            profile for this read length (right).""",
                                            style={"font-style": "italic"},
                                        ),
                                    ],
                                    className="box",
                                    style={"width": "60%"},
                                ),
                                html.Div(
                                    [
                                        html.Img(
                                            id="read_length_metagene_bars_fig",
                                            src="",
                                            style={"margin": "10px"},
                                        ),  # add hoc placement
                                        dcc.Markdown(
                                            id="read_length_metagene_bars_text",
                                            style={"font-size": "large"},
                                        ),
                                    ],
                                    className="box",
                                    style={"width": "40%"},
                                ),
                            ],
                            className="row",
                        ),
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.H2(
                                            """3. Periodicity for all selected samples"""
                                        ),
                                        html.Br(),
                                        dcc.Graph(id="stacked_frames_fig"),
                                    ],
                                    className="box",
                                    style={"width": "40%"},
                                ),
                                # TODO add description ?
                                html.Div(
                                    [
                                        dcc.Input(
                                            id="window_upstream",
                                            type="number",
                                            placeholder="Nucleotides upstream",
                                            min=0,
                                            max=min_upstream,
                                            # value=3,
                                            style={"margin-right": "15px"},
                                        ),
                                        dcc.Input(
                                            id="window_downstream",
                                            type="number",
                                            placeholder="Nucleotides downstream",
                                            min=0,
                                            max=min_downstream,
                                            # value=3,
                                            style={"margin-right": "15px"},
                                        ),
                                        dcc.Input(
                                            id="window_step",
                                            type="number",
                                            placeholder="Step size",
                                            min=1,
                                            max=50,
                                            # value=1,
                                            style={"margin-right": "15px"},
                                        ),
                                        html.Br(),
                                        html.Div(
                                            [
                                                html.Label(
                                                    "Select plot type",
                                                    style={
                                                        "margin-top": "25px",
                                                        "margin-right": "10px",
                                                    },
                                                ),
                                                html.Div(
                                                    [
                                                        dcc.Dropdown(
                                                            id="metagene_type",
                                                            clearable=False,
                                                            searchable=False,
                                                            options=[
                                                                "Discrete",
                                                                "Continuous",
                                                            ],
                                                            value="Discrete",
                                                            style={
                                                                "margin-top": "20px"
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "width": "25%",
                                                        "margin-right": "25px",
                                                    },
                                                ),
                                            ],
                                            className="row",
                                        ),
                                        html.Div(
                                            [html.Img(id="all_metagene_fig", src="")],
                                            id="plot_div",
                                        ),
                                    ],
                                    className="box",
                                    style={"width": "60%"},
                                ),
                            ],
                            className="row",
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


# stacked bars read filtering counts
@app.callback(
    Output("stacked_reads_fig", "figure"),
    [Input("drop_samples", "value"), Input("radio_stacked_reads", "value")],
)
def stacked_reads(selected, zoom):

    selected_samples = selected.split(",")

    stack_cts_order = STACK_CTS_ORDER
    stack_cts_name = STACK_CTS_NAME
    if zoom == 1:  # exclude ribosomal/poor quality
        stack_cts_order = stack_cts_order[:4]
        stack_cts_name = stack_cts_name[:4]
    pal = pal_set3[: len(stack_cts_name)]

    # get differential counts - ascending order
    alignment_diff_counts = get_diff_counts(read_filtering_counts[stack_cts_order])
    df = pd.DataFrame(alignment_diff_counts)
    df.columns = stack_cts_name
    df["Sample"] = read_filtering_counts["Sample"].reset_index(drop=True)
    df = df[df.Sample.isin(selected_samples)]

    fig = (
        px.bar(
            df,
            x="Sample",
            y=stack_cts_name,
            color_discrete_sequence=pal,
            labels={"value": "Reads", "variable": "Filter"},
        )
        .update_layout(
            {"plot_bgcolor": "rgba(0,0,0,0)", "paper_bgcolor": "rgba(0,0,0,0)"}
        )
        .update_layout(legend_traceorder="reversed")
    )
    return fig


# funnel chart
@app.callback(
    Output("funnel_fig", "figure"),
    [
        Input("stacked_reads_fig", "hoverData"),
        Input("radio_stacked_reads", "value"),
        Input("drop_samples", "value"),
    ],
)
def funnel_reads(hoverData, zoom, selected):

    stack_cts_order = STACK_CTS_ORDER
    funnel_cts_name = FUNNEL_CTS_NAME
    if zoom == 1:  # exclude ribosomal/poor quality
        stack_cts_order = stack_cts_order[:4]
        funnel_cts_name = funnel_cts_name[:4]
    pal = pal_set3[: len(stack_cts_order)]

    if hoverData is not None:
        sample_name = hoverData["points"][0]["label"]
        title = f"Filtering steps for {sample_name}"
    else:
        sample_name = selected.split(",")[0]  # pick first - TODO: ordering/grouping
        title = f"Filtering steps for {sample_name}"

    df = read_filtering_counts[read_filtering_counts.Sample == sample_name]
    df = df[stack_cts_order[::-1]]
    df.columns = funnel_cts_name

    fig = (
        go.Figure(
            go.Funnel(
                y=df.columns,
                x=df.values[0],
                textinfo="label+value",
                hoverinfo="text+percent previous+percent initial",
                marker={"color": pal_set3[: len(stack_cts_order)][::-1]},
            )
        )
        .update_yaxes(showticklabels=False)
        .update_layout(
            {"plot_bgcolor": "rgba(0,0,0,0)", "paper_bgcolor": "rgba(0,0,0,0)"}
        )
        .update_layout(
            title_text=title,
            title_x=0.5,
            autosize=True,
            margin=dict(t=100, b=0, l=5, r=5),
        )
    )
    return fig


# available samples
@app.callback(
    [
        Output("dropdown_available_samples", "options"),
        Output("dropdown_available_samples", "value"),
    ],
    Input("drop_samples", "value"),
)
def update_available_samples(selected):
    options = [dict(label=val, value=val) for val in selected.split(",")]
    return (options, options[0]["value"])


# TODO: can we share between callbacks?
# read length distribution min max for sample
@app.callback(
    [
        Output("read_length_slider", "min"),
        Output("read_length_slider", "max"),
        Output("read_length_slider", "value"),
    ],
    [Input("dropdown_available_samples", "value"), Input("drop_unique", "value")],
)
def read_length_min_max(sample, unique):

    unique = ast.literal_eval(unique)
    m_unique = read_length_distributions.is_unique == unique
    m_sample = read_length_distributions.Sample == sample
    df = read_length_distributions[m_sample & m_unique]
    min_length = df.length.min()
    max_length = df.length.max()

    return (min_length, max_length, [min_length, max_length])


# read length distribution
@app.callback(
    Output("length_bars_fig", "figure"),
    [
        Input("dropdown_available_samples", "value"),
        Input("drop_unique", "value"),
        Input("toggle_log", "value"),
        Input("read_length_slider", "value"),
    ],
)
def length_bars(sample, unique, log_y, slider):

    unique = ast.literal_eval(unique)
    pal = pal_bars[2]
    if unique:
        pal = pal_bars[1]
    color_discrete_map = {"Non-periodic": pal, "Periodic": pal_bars[0]}

    m_unique = read_length_distributions.is_unique == unique
    m_sample = read_length_distributions.Sample == sample
    m_min = read_length_distributions.length >= slider[0]
    m_max = read_length_distributions.length <= slider[1]
    m_view = m_unique & m_sample & m_min & m_max
    df = read_length_distributions[m_view].copy()  # WHY
    df["Status"] = "Non-periodic"
    m_sample = lengths_and_offsets.Sample == sample
    m_periodic = lengths_and_offsets.status == "Used for analysis"
    periodic_lengths = lengths_and_offsets[m_sample & m_periodic].length.values
    df.loc[df.length.isin(periodic_lengths), "Status"] = "Periodic"

    fig = (
        px.bar(
            df,
            x="length",
            y="count",
            log_y=log_y,
            color="Status",
            color_discrete_map=color_discrete_map,
            text="length",
            labels={"length": "Length", "count": "Read count"},
        )
        .update_layout(
            {"plot_bgcolor": "rgba(0,0,0,0)", "paper_bgcolor": "rgba(0,0,0,0)"}
        )
        .update_xaxes(showticklabels=False)
        .update_traces(textposition="outside", cliponaxis=False)
    )

    return fig


# read length periodicity
@app.callback(
    [
        Output("read_length_metagene_bars_fig", "src"),
        Output("read_length_metagene_bars_text", "children"),
    ],
    [
        Input("length_bars_fig", "hoverData"),
        Input("dropdown_available_samples", "value"),
    ],
)
def read_length_metagene_bars(hoverData, sample):

    # TODO: how to get a default value, i.e. not hoverData['points'][0] is None?

    m_sample = lengths_and_offsets.Sample == sample
    sample_lengths_and_offsets = lengths_and_offsets[m_sample]

    if hoverData is None:
        length = sample_lengths_and_offsets.length.values.min()
    else:
        length = hoverData["points"][0]["label"]

    metagene_profiles = pd.read_csv(
        filenames.get_metagene_profiles(
            config["riboseq_data"],
            reverse_sample_name_map[sample],
            is_unique=is_unique,
            note=config_note,
        )
    )

    m_length = metagene_profiles["length"] == length
    metagene_profile = metagene_profiles[m_length]
    start_counts, end_counts = get_window_counts(
        metagene_profile,
        0,  # no offset
        read_length_start_upstream,
        read_length_start_downstream,
        read_length_end_upstream,
        read_length_end_downstream,
    )

    # for some read lengths there are no profile...
    # for now, we fill with zeros...
    start_window_size = read_length_start_downstream - read_length_start_upstream + 1
    end_window_size = read_length_end_downstream - read_length_end_upstream + 1
    if not start_counts.any():
        start_counts = np.zeros(start_window_size)
    if not end_counts.any():
        end_counts = np.zeros(end_window_size)

    # otherwise identical low and high ylims makes transformation singular
    ymax = max(max(start_counts), max(end_counts))
    y_max = None
    if ymax == 0:
        ymax = 1
        y_max = ymax
    ymax = 1.1 * ymax

    figsize = (5.5, 4)
    y_label = "Ribo-seq reads"
    fig, axes = metagene_plot_bars(
        start_counts,
        end_counts,
        read_length_start_upstream,
        read_length_start_downstream,
        read_length_end_upstream,
        read_length_end_downstream,
        read_length_step,
        figsize,
        y_label,
        y_max=y_max,
        layout="constrained",
    )

    # selected length and offset if available...
    try:
        m_length = sample_lengths_and_offsets.length == length
        sample_lengths_and_offsets = sample_lengths_and_offsets[m_length]
        offset = sample_lengths_and_offsets.offset.values[0]
        status = sample_lengths_and_offsets.status.values[0]
    except:
        offset = "NA"
        status = "NA"

    # start
    translation_start_position = -1 * read_length_start_upstream
    axes[0].plot(
        [translation_start_position, translation_start_position],
        [0, ymax],
        color="k",
        alpha=0.5,
        linestyle="--",
    )
    # p-site
    if offset != "NA":
        psite = translation_start_position + offset

        axes[0].plot(
            [psite, psite],
            [0, ymax],
            color="k",
            alpha=0.25,
            linestyle="--",
        )
    # end
    translation_end_position = -1 * read_length_end_upstream
    axes[1].plot(
        [translation_end_position, translation_end_position],
        [0, ymax],
        color="k",
        alpha=0.5,
        linestyle="--",
    )

    # remove ticks
    labels = [
        "" if item.get_text() != "0" else "TIS" for item in axes[0].get_xticklabels()
    ]
    axes[0].set_xticklabels(labels)
    labels = [
        "" if item.get_text() != "0" else "TTS" for item in axes[1].get_xticklabels()
    ]
    axes[1].set_xticklabels(labels)

    default_xlabel_start = "Position of P-site relative to start (nt)"
    default_xlabel_end = "Position of P-site relative to stop (nt)"
    axes[0].set_xlabel(default_xlabel_start, fontsize=8)
    axes[1].set_xlabel(default_xlabel_end, fontsize=8)

    # set y ticks
    # axes[0].yaxis.set_major_formatter(ticker.EngFormatter())

    text = (
        f"Read length periodicity\n\n"
        f"---\n\n"
        f"**Selected sample:** {sample}\n\n"
        f"**Length:** {length}\n\n"
        f"**P-site offset:** {offset}\n\n"
        f"**Status:** {status}\n\n"
    )
    out_url = fig_to_uri(fig)

    return (out_url, text)


# percentage of reads in each frame
@app.callback(Output("stacked_frames_fig", "figure"), Input("drop_samples", "value"))
def stacked_frames(selected):

    selected_samples = selected.split(",")
    df = frame_counts[frame_counts.Sample.isin(selected_samples)]

    fig = (
        px.bar(
            df,
            x="Sample",
            y="p_site",
            color_discrete_sequence=pal_frames,
            color="Frame",
            labels={"p_site": "% P-sites (ORF profiles)"},
        )
        .update_layout(
            {"plot_bgcolor": "rgba(0,0,0,0)", "paper_bgcolor": "rgba(0,0,0,0)"}
        )
        .update_layout(legend_traceorder="reversed")
    )

    return fig


# across all samples metagene profile (bars)
@app.callback(
    Output("all_metagene_fig", "src"),
    [
        Input("window_upstream", "value"),
        Input("window_downstream", "value"),
        Input("window_step", "value"),
        Input("metagene_type", "value"),
        Input("drop_samples", "value"),
    ],
)
def all_metagenes(window_upstream, window_downstream, step, plot, selected):

    # instead of defining a default value and override placeholder...
    if window_upstream is None:
        window_upstream = 3
    if window_downstream is None:
        window_downstream = 3
    if step is None:
        step = 1

    start_upstream_window = -window_upstream
    start_downstream_window = window_downstream
    end_upstream_window = -window_downstream
    end_downstream_window = window_upstream

    selected_samples = selected.split(",")

    if plot == "Discrete":
        metagene_profile_start, metagene_profile_end = get_profiles_bars(
            selected_samples,
            start_upstream_window,
            start_downstream_window,
            end_upstream_window,
            end_downstream_window,
        )

        figsize = (8, 6)
        y_label = "Periodic Ribo-seq reads"
        fig, axes = metagene_plot_bars(
            metagene_profile_start,
            metagene_profile_end,
            start_upstream_window,
            start_downstream_window,
            end_upstream_window,
            end_downstream_window,
            step,
            figsize,
            y_label,
        )

        label = "P-site shifted position relative to annotated CDSs (nt)"
        fig.text(0.5, 0.01, label, ha="center")  # adjust placement...
        out_url = fig_to_uri(fig)
    else:  # Continuous
        metagene_profile_start, metagene_profile_end = get_profiles_lines(
            selected_samples,
            start_upstream_window,
            start_downstream_window,
            end_upstream_window,
            end_downstream_window,
        )

        figsize = (8, 6)
        y_label = "Normalised density (AU)"
        fig, axes = metagene_plot_lines(
            metagene_profile_start,
            metagene_profile_end,
            start_upstream_window,
            start_downstream_window,
            end_upstream_window,
            end_downstream_window,
            step,
            figsize,
            y_label,
        )

        label = "P-site shifted position relative to annotated CDSs (nt)"
        fig.text(0.5, 0.01, label, ha="center")  # adjust placement...
        out_url = fig_to_uri(fig)

    return out_url


def main():
    if "DISPLAY" in os.environ and os.environ["DISPLAY"]:
        threading.Timer(
            1, lambda: webbrowser.open_new(f"http://{host}:{port}/")
        ).start()
    app.run(debug=debug, host=host, port=port)


if __name__ == "__main__":
    main()
