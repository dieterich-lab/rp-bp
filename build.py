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
    "pysam",
    "pyfasta",
    "pystan",
    "pybedtools",
    "pyyaml",
    "-e git+https://bitbucket.org/bmmalone/misc.git#egg=misc[bio]",
    "-e ."
]

clean_requirements = [
    "cython",
    "numpy",
    "scipy",
    "pandas",
    "joblib",
    "docopt",
    "tqdm",
    "pysam",
    "pyfasta",
    "pystan",
    "pybedtools",
    "pyyaml",
    "misc[bio]",
    "rpbp"
]

stan_model_files = [
    os.path.join(os.getcwd(), "models", "metagene-periodicity", "no-periodicity.stan"),
    os.path.join(os.getcwd(), "models", "metagene-periodicity", "start-high-high-low.stan"),
    os.path.join(os.getcwd(), "models", "metagene-periodicity", "start-high-low-high.stan"),
    os.path.join(os.getcwd(), "models", "metagene-periodicity", "start-high-low-low.stan"),
    os.path.join(os.getcwd(), "models", "orf-translation", "gaussian-naive-bayes.stan"),
    os.path.join(os.getcwd(), "models", "orf-translation", "periodic-cauchy-mixture.stan"),
    os.path.join(os.getcwd(), "models", "orf-translation", "zero-inflated-periodic-cauchy-mixture.stan")
]

stan_pickle_files = [
    os.path.join(os.getcwd(), "models", "metagene-periodicity", "no-periodicity.pkl"),
    os.path.join(os.getcwd(), "models", "metagene-periodicity", "start-high-high-low.pkl"),
    os.path.join(os.getcwd(), "models", "metagene-periodicity", "start-high-low-high.pkl"),
    os.path.join(os.getcwd(), "models", "metagene-periodicity", "start-high-low-low.pkl"),
    os.path.join(os.getcwd(), "models", "orf-translation", "gaussian-naive-bayes.pkl"),
    os.path.join(os.getcwd(), "models", "orf-translation", "periodic-cauchy-mixture.pkl"),
    os.path.join(os.getcwd(), "models", "orf-translation", "zero-inflated-periodic-cauchy-mixture.pkl")
]

def install(openblas_install_path):
    global requirements

    # first, check to see if OpenBLAS is already installed
    open_blas_symlink_path = os.path.join(openblas_install_path, 'libopenblas.so')
    open_blas_0_symlink_path = os.path.join(openblas_install_path, 'libopenblas.so.0')

    if not os.path.exists(open_blas_symlink_path):
        msg = "The OpenBLAS libraries were not found. They will be downloaded and compiled."
        logging.info(msg)

        # then download and make OpenBLAS
        cmd = "git clone https://github.com/xianyi/OpenBLAS"
        subprocess.call(cmd, shell=True)

        cmd = "cd OpenBLAS && make FC=gfortran"
        subprocess.call(cmd, shell=True)

        # and create symlinks to the read files to the desired library path
        open_blas_library_path = os.path.join(os.getcwd(), 'OpenBLAS', "libopenblas.so")
        open_blas_0_library_path = os.path.join(os.getcwd(), 'OpenBLAS', "libopenblas.so.0")

        msg = "Symlinks will be created to the OpenBLAS libraries."
        logging.info(msg)

        if os.path.exists(open_blas_symlink_path):
            os.remove(open_blas_symlink_path)
        os.symlink(open_blas_library_path, open_blas_symlink_path)

        if os.path.exists(open_blas_0_symlink_path):
            os.remove(open_blas_0_symlink_path)
        os.symlink(open_blas_0_library_path, open_blas_0_symlink_path)

        

    else:
        msg = "The OpenBLAS libraries were found and will be used for linking."
        logging.info(msg)

    # additionally, check if we are working in a virtual environment
    virtual_env_path = os.getenv('VIRTUAL_ENV')
    if virtual_env_path is not None:
        msg = "Symlinks will also be copied to the virtual environment lib directory"
        logging.info(msg)

        virtual_env_lib_path = os.path.join(virtual_env_path, 'lib', 'libopenblas.so')
        virtual_env_0_lib_path = os.path.join(virtual_env_path, 'lib', 'libopenblas.so.0')
        
        if os.path.exists(virtual_env_lib_path):
            os.remove(virtual_env_lib_path)
        os.symlink(open_blas_library_path, virtual_env_lib_path)

        if os.path.exists(virtual_env_0_lib_path):
            os.remove(virtual_env_0_lib_path)
        os.symlink(open_blas_0_library_path, virtual_env_0_lib_path)


    # either way, make it clear they need to be in a place where the linker can find them
    msg = ("The OpenBLAS libraries will be used for the python libraries. Please ensure "
        "the following path is in the LD_LIBRARY_PATH variable (or analog for your "
        "operating system).\n\n{}\n\nFor bash, this line can be added to the ~/.bashrc "
        "file to ensure correct compilation:\n\nexport LD_LIBRARY_PATH={}:$LD_LIBRARY_PATH\n\n".format(
        openblas_install_path, openblas_install_path))
    logging.info(msg)

    
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

    # remove OpenBLAS
    if os.path.exists('OpenBLAS'):
        msg = "Removing the local OpenBLAS installation"
        logging.info(msg)

        shutil.rmtree('OpenBLAS')

        msg = "Removing OpenBLAS symlinks"
        logging.info(msg)

        open_blas_symlink_path = os.path.join(openblas_install_path, 'libopenblas.so')
        open_blas_0_symlink_path = os.path.join(openblas_install_path, 'libopenblas.so.0')

        if os.path.lexists(open_blas_symlink_path):
            os.unlink(open_blas_symlink_path)
        
        if os.path.lexists(open_blas_0_symlink_path):
            os.unlink(open_blas_0_symlink_path)

        
        # additionally, check if we are working in a virtual environment
        virtual_env_path = os.getenv('VIRTUAL_ENV')
        if virtual_env_path is not None:
            msg = "Symlinks will also be removed from the virtual environment lib directory"
            logging.info(msg)
        
            virtual_env_lib_path = os.path.join(virtual_env_path, 'lib', 'libopenblas.so')
            virtual_env_0_lib_path = os.path.join(virtual_env_path, 'lib', 'libopenblas.so.0')
            
            if os.path.lexists(virtual_env_lib_path):
                os.unlink(virtual_env_lib_path)

            if os.path.lexists(virtual_env_0_lib_path):
                os.unlink(virtual_env_0_lib_path)
    else:
        msg = ("The OpenBLAS libraries appear to have not been installed by this package. "
            "They will not be removed.")
        logging.info(msg)

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
