
# Installing the Rp-Bp pipeline

This document describes detailed installation instructions for the Rp-Bp pipeline. These steps have been primarily tested on Ubuntu.

<a id="toc"></a>

* [Prerequisites](#prerequisites)
* [Installation](#installation)
    * [Virtual environment installation](#virtual-environment-installation)
    * [Anaconda installation](#anaconda-installation)
    * [Local installation](#local-installation)
* [Uninstallation](#uninstallation)

---

<a id="prerequisites"></a>

## Prerequisites

The pipelines make use of a number of standard bioinformatics tools. All of these must be installed and available on the `$PATH` for the pipeline to work correctly. All of the pipeline scripts check that the required programs are available before executing. If any required programs cannot be found, the script prints an error message about the missing program and does not continue. The versions used during development and testing are specified below. For most of the tools, any recent version should be sufficient. If problems arise, though, please use the version indicated below

* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), version 2.2.6
* [flexbar](https://github.com/seqan/flexbar), version 2.5
* [SAMtools](http://www.htslib.org/), version 1.2
* [STAR](https://github.com/alexdobin/STAR), version 2.4.1d

#### Helper script

While not officially supported, the script
[here](https://gist.github.com/bmmalone/43752eba0af97d1085eef7db033309d0)
downloads and compiles the necessary prerequisites in `$HOME/install` and
installs them in `$HOME/local`.

It uses `apt-get` to install the Intel thread building blocks and CMake. The
direct Rp-Bp dependencies mentioned above are simply downloaded and placed in
the appropriate location. The relevant location (`$HOME/local/bin`) must be
added to the `$PATH` before the executable files can be found. `sudo` permissions
are required for the calls to `apt-get`.

#### OpenBLAS

If installation fails due to missing `OpenBLAS` dependencies for `scipy`, please follow the instructions [here](https://gist.github.com/bmmalone/1b5f9ff72754c7d4b313c0b044c42684).

[Back to top](#toc)

<a id='installation'></a>

## Installation

The installation currently uses the requirement file `requirements.txt` to install a number of required packages, as explained below. 

<a id='virtual-environment-installation'></a>

## Virtual environment installation

We encourage installing Rp-Bp and dependencies in a virtual environment. The `venv` module provides support for creating environments with their own site directories. See [venv](https://docs.python.org/3/library/venv.html) for more information about Python virtual environments. To create a virtual environment:

```
python3 -m venv /path/to/virtual/environment
```
To activate the new virtual environment and install Rp-Bp:

```
# Activate the new virtual environment.
source /path/to/virtual/environment/bin/activate

# If necessary, upgrade pip and additional packages (such as setuptools if installing in editable mode).
pip install --upgrade pip

# Clone the git repository
git clone git@github.com:dieterich-lab/rp-bp.git

# Change to the rp-bp directory and build the package. 
cd rp-bp && pip install --verbose -r requirements.txt

```

The build process includes compiling several libraries for optimized numerical calculations. Due to the optimized nature of these libraries, the initial installation can take up to an hour.

Due to the `--verbose` flag, much debugging information is printed. In some cases, building some libraries may initially fail; this is due to a known issue with dependency-handling in [pip](https://pip.pypa.io/en/stable/reference/pip_install/#installation-order). The following output (or something similar depending on packages already available) indicates that installation succeeded:

```
Successfully installed cython-0.25.1 docopt-0.6.2 joblib-0.10.3 numpy-1.11.2 pandas-0.19.0 patsy-0.4.1 psutil-4.4.2 pyfasta-0.5.2 pysam-0.9.1.4 pystan-2.12.0.0 python-dateutil-2.5.3 pytz-2016.7 pyyaml-3.12 rpbp-1.0 scipy-0.18.1 six-1.10.0 statsmodels-0.6.1 tqdm-4.9.0
Cleaning up...
```

To use the programs in the future, use the `source` command to activate the virtual environment, as described above.

[Back to top](#toc)

<a id='anaconda-installation'></a>

## Anaconda installation

The package can also be installed within an [anaconda](https://www.continuum.io/) environment. 

   * The version of gcc included by default may not properly compile the Stan model files, so it must be updated first.
   * In addition, `llvm` may have to be installed (required for `pymisc`):
    
```
# Update gcc if necessary.
conda install -c salford_systems libgcc-6=6.2.0

# Install llvm
conda install -c anaconda llvm
conda install -c numba llvmlite

# Create the anaconda environment.
conda create -n my_new_environment python=3.5 anaconda

# Activate the new environment.
source activate my_new_environment

# Clone the git repository
git clone git@github.com:dieterich-lab/rp-bp.git

# Change to the rp-bp directory and build the package. 
cd rp-bp && pip install --verbose -r requirements.txt
```

[Back to top](#toc)

<a id='local-installation'></a>

## Local installation

We recommend installing Rp-Bp in a virtual environment as described [above](#virtual-environment-installation). However if needed, the following instructions can be used to install the package without sudo access in a user's home directory.

The commands below are presumably executed in a directory like `$HOME/install`. They install the python executables into `$HOME/local/bin` (or wherever the `prefix` option is located), so that directory must be in the `$PATH`.
This can be accomplished by adding a line like `export PATH=$HOME/local/bin:$PATH` in the file `.bashrc` on Ubuntu.

   * It is very important that pip and wheel are upgraded!

   * The `--user` option can also be given to `pip3` for installation within the user's home directory. The installation is then compatible with system-wide Python3 installations without requiring sudo access.

```
# Lines to add to .bashrc for the installation process.
export PATH=$HOME/local/bin:$PATH

### Downloading and installing the required software.

# Download, extract and install (locally) Python 3.  
wget https://www.python.org/ftp/python/3.5.1/Python-3.5.1.tgz && tar -xvf Python-3.5.1.tgz && cd Python-3.5.1 && ./configure --prefix=$HOME/local --with-ensurepip=upgrade && make && make install && cd ..

# Upgrade pip. 
pip3 install --upgrade pip wheel

# Clone the git repository. 
git clone git@github.com:dieterich-lab/rp-bp.git
    
# Change into the rp-bp directory and build the package... 
cd rp-bp && pip3 install --verbose -r requirements.txt

# ...or install the package in the user's home directory
cd rp-bp && pip3 install --verbose --user  -r requirements.txt
```

The build process includes compiling several libraries for optimized numerical calculations. Due to the optimized nature of these libraries, the initial installation can take up to an hour.

Due to the `--verbose` flag, much debugging information is printed. In some cases, building some libraries may initially fail; this is due to a known issue with dependency-handling in [pip](https://pip.pypa.io/en/stable/reference/pip_install/#installation-order). The following output (or something similar depending on packages already available) indicates installation succeeded:

```
Successfully installed cython-0.25.1 docopt-0.6.2 joblib-0.10.3 numpy-1.11.2 pandas-0.19.0 patsy-0.4.1 psutil-4.4.2 pyfasta-0.5.2 pysam-0.9.1.4 pystan-2.12.0.0 python-dateutil-2.5.3 pytz-2016.7 pyyaml-3.12 rpbp-1.0 scipy-0.18.1 six-1.10.0 statsmodels-0.6.1 tqdm-4.9.0
Cleaning up...
```

[Back to top](#toc)

<a id='uninstallation'></a>

## Uninstallation

The `rpbp` python and required packages can be removed with pip:

``pip3 uninstall misc bio-utils riboutils rpbp``

If the packages were installed in a dedicated virtual environment, this environment can simply be removed.