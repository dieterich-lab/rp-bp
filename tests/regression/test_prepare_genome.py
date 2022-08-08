
import pytest
import gzip


def get_lines(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, 'rb') as f:
            return f.readlines()
    with open(filename, 'r') as f:
        return f.readlines()


def files_match(file, ref_file):
    try:
        assert get_lines(file) == get_lines(ref_file)
    except AssertionError as e:
        e.args += ("File:", file, "Reference file:", ref_file)
        raise
    return True


def test_output(get_files):
    for file, ref_file in get_files:
        assert files_match(file, ref_file)
