"""
    Regression test using the c-elegans reference dataset.
"""

import logging
import pytest
import gzip
import pandas as pd

logger = logging.getLogger(__name__)


def to_df(filename):
    args = {}
    if filename.endswith(".bed.gz"):
        args = {"sep": "\t"}
    elif filename.endswith(".mtx.gz"):
        args = {"sep": " ", "comment": "%", "header": None}
    return pd.read_csv(filename, **args)


# test output of `prepare-rpbp-genome`
def test_pipeline_part1(getf_genome):
    for file, ref_file in getf_genome:
        msg = f"Comparing {file} and {ref_file}"
        logger.info(msg)
        pd.testing.assert_frame_equal(to_df(file), to_df(ref_file))


# test output of `run-all-rpbp-instances`
@pytest.mark.depends(on=["test_pipeline_part1"])
def test_pipeline_part2(getf_pipeline):
    for file, ref_file in getf_pipeline:
        msg = f"Comparing {file} and {ref_file}"
        logger.info(msg)
        pd.testing.assert_frame_equal(to_df(file), to_df(ref_file))
