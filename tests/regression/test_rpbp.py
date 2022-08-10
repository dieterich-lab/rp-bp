"""
    Regression test using the c-elegans reference dataset.
"""

import logging
import pytest
import gzip

logger = logging.getLogger(__name__)


def get_lines(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, "rb") as f:
            return f.readlines()
    with open(filename, "r") as f:
        return f.readlines()


def files_match(file, ref_file):
    try:
        assert get_lines(file) == get_lines(ref_file)
    except AssertionError as e:
        e.args += ("File:", file, "Reference file:", ref_file)
        raise
    return True


# test output of `prepare-rpbp-genome`
def test_pipeline_part1(getf_genome):
    for file, ref_file in getf_genome:
        msg = f"Comparing {file} and {ref_file}"
        logger.info(msg)
        assert files_match(file, ref_file)
        

# test output of `run-all-rpbp-instances`
def test_pipeline_part2(getf_pipeline):
    for file, ref_file in getf_pipeline:
        msg = f"Comparing {file} and {ref_file}"
        logger.info(msg)
        assert files_match(file, ref_file)
