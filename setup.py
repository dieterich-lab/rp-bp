#! /usr/bin/env python3

import os
import logging
import importlib

from pathlib import Path
from cmdstanpy import CmdStanModel

from setuptools import setup
from setuptools.command.install import install as install
from setuptools.command.develop import develop as develop

logger = logging.getLogger(__name__)

# defaults.py
cpp_options = {"STAN_THREADS": "TRUE"}

stan_model_files = [
    Path("nonperiodic", "no-periodicity.stan"),
    Path("nonperiodic", "start-high-high-low.stan"),
    Path("nonperiodic", "start-high-low-high.stan"),
    Path("periodic", "start-high-low-low.stan"),
    Path("untranslated", "gaussian-naive-bayes.stan"),
    Path("translated", "periodic-gaussian-mixture.stan"),
]


def _post_install(recompile):

    import site

    importlib.reload(site)

    import rpbp.ribo_utils.filenames as filenames
    import pbiotools.misc.shell_utils as shell_utils

    # instantiate models under "models_base"
    models_base = filenames.get_default_models_base()
    src_smfs = [Path("rpbp_models", s) for s in stan_model_files]
    dest_smfs = [Path(models_base, s) for s in stan_model_files]
    for src_smf, dest_smf in zip(src_smfs, dest_smfs):
        dest_smf.parent.mkdir(parents=True, exist_ok=True)
        dest_smf.write_text(src_smf.read_text())
        CmdStanModel(stan_file=dest_smf, compile=recompile, cpp_options=cpp_options)

    # check dependencies
    programs = ["flexbar"]
    shell_utils.check_programs_exist(
        programs, raise_on_error=False, package_name="flexbar", logger=logger
    )

    programs = ["STAR"]
    shell_utils.check_programs_exist(
        programs, raise_on_error=False, package_name="STAR", logger=logger
    )

    programs = ["bowtie2", "bowtie2-build-s"]
    shell_utils.check_programs_exist(
        programs, raise_on_error=False, package_name="bowtie2", logger=logger
    )

    programs = ["samtools"]
    shell_utils.check_programs_exist(
        programs, raise_on_error=False, package_name="SAMtools", logger=logger
    )


class SetupInstall(install):

    user_options = install.user_options + [
        ("force-recompile", None, "Set this flag to recompile the Stan models"),
    ]

    def initialize_options(self):
        install.initialize_options(self)
        self.force_recompile = None

    def finalize_options(self):
        install.finalize_options(self)

    def run(self):
        recompile = True
        if self.force_recompile:
            recompile = "force"

        level = logging.getLevelName("INFO")
        logging.basicConfig(level=level, format="%(levelname)-8s : %(message)s")

        install.run(self)
        # skip if RTD
        if not os.environ.get("READTHEDOCS") == "True":
            _post_install(recompile)


class SetupDevelop(develop):

    user_options = develop.user_options + [
        ("force-recompile", None, "Set this flag to recompile the Stan models"),
    ]

    def initialize_options(self):
        develop.initialize_options(self)
        self.force_recompile = None

    def finalize_options(self):
        develop.finalize_options(self)

    def run(self):
        recompile = True
        if self.force_recompile:
            recompile = "force"

        level = logging.getLevelName("INFO")
        logging.basicConfig(level=level, format="%(levelname)-8s : %(message)s")

        develop.run(self)
        # skip if RTD
        if not os.environ.get("READTHEDOCS") == "True":
            _post_install(recompile)


setup(cmdclass={"install": SetupInstall, "develop": SetupDevelop})
