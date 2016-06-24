#! /usr/bin/env python3

import argparse
import logging
import shutil
import subprocess

import os

default_open_blas_install_dir = os.path.join(os.path.expanduser('~'), 'local', 'lib')

# The list of requirements for this package. They will be installed,
# one at a time, by pip3 in the order they appear in the list.
install_requirements = [
    "cython",
    "numpy",
    "scipy",
    "pandas",
    "joblib",
    "docopt",
    "tqdm",
    "statsmodels",
    "pysam",
    "pyfasta",
    "pystan",
    "pybedtools",
    "pyyaml",
    "git+https://bitbucket.org/bmmalone/misc.git#egg=misc[bio]",
    "git+https://github.com/dieterich-lab/riboseq-utils.git#egg=riboutils",
    "-e .[analysis]"
]

clean_requirements = [
    "cython",
    "numpy",
    "scipy",
    "pandas",
    "joblib",
    "docopt",
    "tqdm",
    "statsmodels",
    "pysam",
    "pyfasta",
    "pystan",
    "pybedtools",
    "pyyaml",
    "misc[bio]",
    "riboutils",
    "rpbp[analysis]"
]

stan_model_files = [
    os.path.join(os.getcwd(), "models", "nonperiodic", "no-periodicity.stan"),
    os.path.join(os.getcwd(), "models", "nonperiodic", "start-high-high-low.stan"),
    os.path.join(os.getcwd(), "models", "nonperiodic", "start-high-low-high.stan"),
    os.path.join(os.getcwd(), "models", "periodic", "start-high-low-low.stan"),
    os.path.join(os.getcwd(), "models", "untranslated", "gaussian-naive-bayes.stan"),
    os.path.join(os.getcwd(), "models", "translated", "periodic-cauchy-mixture.stan"),
    os.path.join(os.getcwd(), "models", "translated", "zero-inflated-periodic-cauchy-mixture.stan")
]

stan_pickle_files = [
    os.path.join(os.getcwd(), "models", "nonperiodic", "no-periodicity.pkl"),
    os.path.join(os.getcwd(), "models", "nonperiodic", "start-high-high-low.pkl"),
    os.path.join(os.getcwd(), "models", "nonperiodic", "start-high-low-high.pkl"),
    os.path.join(os.getcwd(), "models", "periodic", "start-high-low-low.pkl"),
    os.path.join(os.getcwd(), "models", "untranslated", "gaussian-naive-bayes.pkl"),
    os.path.join(os.getcwd(), "models", "translated", "periodic-cauchy-mixture.pkl"),
    os.path.join(os.getcwd(), "models", "translated", "zero-inflated-periodic-cauchy-mixture.pkl")
]

def install(openblas_install_path):
    global requirements
    
    option = "install --no-binary :all:"
    for r in install_requirements:
        cmd = "pip3 {} {}".format(option, r)
        subprocess.call(cmd, shell=True)

    # and now pickle the stan models
    for stan, pickle in zip(stan_model_files, stan_pickle_files):
        cmd = "pickle-stan {} {}".format(stan, pickle)
        subprocess.call(cmd, shell=True)


def clean(openblas_install_path):
    global requirements

    option = "uninstall --yes"
    for r in clean_requirements:
        cmd = "pip3 {} {}".format(option, r)
        subprocess.call(cmd, shell=True)

    # and remove the pickle files
    for pickle in stan_pickle_files:
        if os.path.exists(pickle):
            os.remove(pickle)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script is used to install or uninstall the packages required "
        "for the Rp-Bp and Rp-chi ribosome profiling prediction pipelines.")
    parser.add_argument('--clean', help="If this flag is given, then the packages will "
        "be uninstalled rather than installed.", action='store_true')
    parser.add_argument('--openblas-install-path', help="The location for the symlinks for "
        "the OpenBLAS library.\n\n** This path MUST be in the LD_LIBRARY_PATH for correct "
        "installation. **", default=default_open_blas_install_dir)
    args = parser.parse_args()

    level = logging.getLevelName("INFO")
    logging.basicConfig(level=level,
            format='%(levelname)-8s : %(message)s')

    if args.clean:
        clean(args.openblas_install_path)
    else:
        install(args.openblas_install_path)

if __name__ == '__main__':
    main()
