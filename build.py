#! /usr/bin/env python3

import argparse
import logging
import shutil
import subprocess

import os

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
    "-e git+https://bitbucket.org/bmmalone/misc.git#egg=misc[bio]", 
        # the "-e" seems to be necessary to grab subfolders. I do not
        # understand this, but it seems to work
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
    os.path.join(os.getcwd(), "models", "translated", "periodic-gaussian-mixture.stan")
    #os.path.join(os.getcwd(), "models", "translated", "periodic-cauchy-mixture.stan"),
    #os.path.join(os.getcwd(), "models", "translated", "zero-inflated-periodic-cauchy-mixture.stan")
]

stan_pickle_files = [
    os.path.join(os.getcwd(), "models", "nonperiodic", "no-periodicity.pkl"),
    os.path.join(os.getcwd(), "models", "nonperiodic", "start-high-high-low.pkl"),
    os.path.join(os.getcwd(), "models", "nonperiodic", "start-high-low-high.pkl"),
    os.path.join(os.getcwd(), "models", "periodic", "start-high-low-low.pkl"),
    os.path.join(os.getcwd(), "models", "untranslated", "gaussian-naive-bayes.pkl"),
    os.path.join(os.getcwd(), "models", "translated", "periodic-gaussian-mixture.pkl")
    #os.path.join(os.getcwd(), "models", "translated", "periodic-cauchy-mixture.pkl"),
    #os.path.join(os.getcwd(), "models", "translated", "zero-inflated-periodic-cauchy-mixture.pkl")
]

def check_programs_exist(programs, package_name):
    """ This function checks that all of the programs in the list cam be
        called from python. After checking all of the programs, a message 
        is printed saying which programs could not be found and the package
        where they should be located.

        Internally, this program uses shutil.which, so see the documentation
        for more information about the semantics of calling.

        Arguments:
            programs (list of string): a list of programs to check

        Returns:
            None
    """

    missing_programs = []
    for program in programs:
        exe_path = shutil.which(program)

        if exe_path is None:
            missing_programs.append(program)

    if len(missing_programs) > 0:
        missing_programs_str = ' '.join(missing_programs)
        msg = "Could not find the following programs: {}".format(missing_programs_str)
        print(msg)

        msg = ("Please ensure the {} package is installed before using the Rp-Bp "
            "pipeline.".format(package_name))
        print(msg)


def install():
    global requirements
    
    option = "install --no-binary :all:"
    for r in install_requirements:
        cmd = "pip3 {} {}".format(option, r)
        subprocess.call(cmd, shell=True)

    # and now pickle the stan models
    for stan, pickle in zip(stan_model_files, stan_pickle_files):
        cmd = "pickle-stan {} {}".format(stan, pickle)
        subprocess.call(cmd, shell=True)

    # finally, check for the prerequisite programs

    programs = ['flexbar']
    check_programs_exist(programs, 'flexbar')

    programs = ['STAR']
    check_programs_exist(programs, 'STAR')

    programs = ['bowtie2', 'bowtie2-build-s']
    check_programs_exist(programs, 'bowtie2')

    programs = ['intersectBed', 'bedToBam', 'fastaFromBed']
    check_programs_exist(programs, 'bedtools')

    programs = ['samtools']
    check_programs_exist(programs, 'SAMtools')

    programs = ['gffread']
    check_programs_exist(programs, 'cufflinks')

def clean():
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
    args = parser.parse_args()

    level = logging.getLevelName("INFO")
    logging.basicConfig(level=level,
            format='%(levelname)-8s : %(message)s')

    if args.clean:
        clean()
    else:
        install()

if __name__ == '__main__':
    main()
