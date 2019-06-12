#! /usr/bin/env python3

import importlib
import logging
import os
import subprocess

from setuptools import setup
from setuptools.command.install import install as install
from setuptools.command.develop import develop as develop


logger = logging.getLogger(__name__)


stan_model_files = [
    os.path.join("nonperiodic", "no-periodicity.stan"),
    os.path.join("nonperiodic", "start-high-high-low.stan"),
    os.path.join("nonperiodic", "start-high-low-high.stan"),
    os.path.join("periodic", "start-high-low-low.stan"),
    os.path.join("untranslated", "gaussian-naive-bayes.stan"),
    os.path.join("translated", "periodic-gaussian-mixture.stan")
]


stan_pickle_files = [
    os.path.join("nonperiodic", "no-periodicity.pkl"),
    os.path.join("nonperiodic", "start-high-high-low.pkl"),
    os.path.join("nonperiodic", "start-high-low-high.pkl"),
    os.path.join("periodic", "start-high-low-low.pkl"),
    os.path.join("untranslated", "gaussian-naive-bayes.pkl"),
    os.path.join("translated", "periodic-gaussian-mixture.pkl")
]


def _pickle_it(stan, pickle):

    import shlex

    dirname = os.path.dirname(pickle)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    cmd = "pickle-stan {} {}".format(shlex.quote(stan), shlex.quote(pickle))
    logging.info(cmd)
    subprocess.call(cmd, shell=True)


def _post_install(force_recompile):

    import site
    importlib.reload(site)

    import pbio.ribo.ribo_filenames as filenames
    import pbio.misc.shell_utils as shell_utils
    
    smf = [os.path.join("rpbp_models", s) for s in stan_model_files]

    models_base = filenames.get_default_models_base()
    spf = [os.path.join(models_base, s) for s in stan_pickle_files]

    # Compile and pickle the Stan models
    if force_recompile:
        for stan, pickle in zip(smf, spf):
            _pickle_it(stan, pickle)
    else:  # default
        for stan, pickle in zip(smf, spf):
            if os.path.exists(pickle):
                msg = "A model already exists at: {}. Skipping.".format(pickle)
                logging.warning(msg)
                continue
            _pickle_it(stan, pickle)

    # Check for the prerequisite programs
    programs = ['flexbar']
    shell_utils.check_programs_exist(programs, raise_on_error=False,
                                     package_name='flexbar', logger=logger)
        
    programs = ['STAR']
    shell_utils.check_programs_exist(programs, raise_on_error=False,
                                     package_name='STAR', logger=logger)

    programs = ['bowtie2', 'bowtie2-build-s']
    shell_utils.check_programs_exist(programs, raise_on_error=False,
                                     package_name='bowtie2', logger=logger)

    programs = ['samtools']
    shell_utils.check_programs_exist(programs, raise_on_error=False,
                                     package_name='SAMtools', logger=logger)


class SetupInstall(install):

    user_options = install.user_options + [
        ('force-recompile', None, 'Set this flag to recompile the Stan models'),
    ]

    def initialize_options(self):
        install.initialize_options(self)
        self.force_recompile = None

    def finalize_options(self):
        install.finalize_options(self)

    def run(self):
        force_recompile = self.force_recompile  # 0 or 1

        level = logging.getLevelName("INFO")
        logging.basicConfig(level=level,
                            format='%(levelname)-8s : %(message)s')

        install.run(self)
        # skip if RTD
        if not os.environ.get('READTHEDOCS') == 'True':
            _post_install(force_recompile)


class SetupDevelop(develop):

    user_options = develop.user_options + [
        ('force-recompile', None, 'Set this flag to recompile the Stan models'),
    ]

    def initialize_options(self):
        develop.initialize_options(self)
        self.force_recompile = None

    def finalize_options(self):
        develop.finalize_options(self)

    def run(self):
        force_recompile = self.force_recompile  # 0 or 1

        level = logging.getLevelName("INFO")
        logging.basicConfig(level=level,
                            format='%(levelname)-8s : %(message)s')

        develop.run(self)
        # skip if RTD
        if not os.environ.get('READTHEDOCS') == 'True':
            _post_install(force_recompile)


setup(
    cmdclass={
        'install': SetupInstall,
        'develop': SetupDevelop
    }
)
