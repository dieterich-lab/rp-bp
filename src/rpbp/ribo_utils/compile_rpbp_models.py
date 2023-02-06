#! /usr/bin/env python3

# Note: not importing pbiotools here to allow the build of the bioconda rp-bp recipe
# to compile the Stan models without needing to install pbiotools

import argparse

from pathlib import Path
from cmdstanpy import CmdStanModel
from rpbp.defaults import model_inst_options

import rpbp.ribo_utils.filenames as filenames

stan_model_files = [
    Path("nonperiodic", "no-periodicity.stan"),
    Path("nonperiodic", "start-high-high-low.stan"),
    Path("nonperiodic", "start-high-low-high.stan"),
    Path("periodic", "start-high-low-low.stan"),
    Path("untranslated", "gaussian-naive-bayes.stan"),
    Path("translated", "periodic-gaussian-mixture.stan"),
]


def _stan_exes_found(models_base) -> bool:
    for stan_model_file in stan_model_files:
        stan_exe_file = Path(models_base, stan_model_file).with_suffix("")
        if not stan_exe_file.is_file():
            return False
    return True


def compile(force_compile: bool = False):
    models_base = filenames.get_default_models_base()

    if _stan_exes_found(models_base) and not force_compile:
        return

    cpp_options = None
    if model_inst_options["stan_threads"]:
        cpp_options = {"STAN_THREADS": "TRUE"}

    for stan_model_file in stan_model_files:
        stan_file = Path(models_base, stan_model_file)
        CmdStanModel(stan_file=stan_file, compile="force", cpp_options=cpp_options)


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Instantiate the Rp-Bp Stan models.
                    The constructor will compile the models if no exe files exist.
                    Default options for C++ compiler are taken from rpbp.defaults.""",
    )
    parser.add_argument(
        "--force",
        help="""Force compilation, even if existing executables are found.""",
        action="store_true",
    )
    args = parser.parse_args()

    compile(args.force)


if __name__ == "__main__":
    main()
