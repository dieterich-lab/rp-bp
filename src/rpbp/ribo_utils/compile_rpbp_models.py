#! /usr/bin/env python3

import logging
import argparse

from pathlib import Path
from cmdstanpy import CmdStanModel
from rpbp.defaults import model_inst_options

import rpbp.ribo_utils.filenames as filenames

import pbiotools.misc.logging_utils as logging_utils

logger = logging.getLogger(__name__)

stan_model_files = [
    Path("nonperiodic", "no-periodicity.stan"),
    Path("nonperiodic", "start-high-high-low.stan"),
    Path("nonperiodic", "start-high-low-high.stan"),
    Path("periodic", "start-high-low-low.stan"),
    Path("untranslated", "gaussian-naive-bayes.stan"),
    Path("translated", "periodic-gaussian-mixture.stan"),
]


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Instantiate the Rp-Bp Stan models. The constructor 
                    will compile the models as needed, i.e. if no exe files 
                    exists or if the Stan files are newer than the exe files.
                    Default options for C++ compiler are taken from rpbp.defaults.""",
    )

    parser.add_argument("--force", help="""Force compilation, even if existing 
                        executables are found.""", action="store_true")

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    recompile = True
    if args.force:
        recompile = "force"
    
    cpp_options = None
    if model_inst_options["stan_threads"]:
        cpp_options = {"STAN_THREADS": "TRUE"}
    
    # instantiate models under "models_base"
    models_base = filenames.get_default_models_base()
    src_smfs = [Path("rpbp_models", s) for s in stan_model_files]
    dest_smfs = [Path(models_base, s) for s in stan_model_files]
    for src_smf, dest_smf in zip(src_smfs, dest_smfs):
        dest_smf.parent.mkdir(parents=True, exist_ok=True)
        if dest_smf.is_file() and not recompile ==  "force":
            msg = f"A model file already exists under {dest_smf}... " \
                  "to force recompilation, use [--force]."
            logger.warning(msg)
            continue
        # copying the Stan model file will anyway trigger recompilation
        dest_smf.write_text(src_smf.read_text())
        CmdStanModel(stan_file=dest_smf, 
                     compile=recompile,
                     cpp_options=cpp_options)
        

if __name__ == "__main__":
    main()
    
