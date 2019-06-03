
# Installing the Rp-Bp pipeline

This document describes detailed installation instructions for the Rp-Bp pipeline. These steps have been primarily tested on Ubuntu.

<a id="toc"></a>

* [Prerequisites](#prerequisites)
* [Installation](#installation)
    * [Virtual environment installation](#virtual-environment-installation)
    * [Anaconda installation](#anaconda-installation)
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

The `rpbp` package can be installed in editable mode by adding the `-e` option to 
the `requirements.txt` file. 

The `--force-recompile` flag can also be passed via `install-option` in the `requirements.txt` file 
to force the recompilation of the Stan models. This flag does not normally have to be set, unless one 
has to re-install `rpbp` (*e.g* after upgrade or changing versions of Pystan). 


<a id='virtual-environment-installation'></a>

## Installation (virtual environment)

To install `rpbp` and dependencies, first create a virtual environment:
 
```
python3 -m venv /path/to/virtual/environment
```

For information about Python virtual environments, see the [venv](https://docs.python.org/3/library/venv.html) documentation.
To activate the new virtual environment and install `rpbp`:

```
# Activate the new virtual environment.
source /path/to/virtual/environment/bin/activate

# If necessary, upgrade pip and wheel or additional packages (such as setuptools if installing in editable mode).
pip install --upgrade pip setuptools wheel

# Clone the git repository
git clone git@github.com:dieterich-lab/rp-bp.git
cd rp-bp

pip --verbose install -r requirements.txt

```

[Back to top](#toc)

## Anaconda installation

The package can also be installed within an [anaconda](https://www.continuum.io/) environment. 

   * The version of gcc included by default may not properly compile the Stan model files, so it must be updated first.
   * In addition, `llvm` and `libyaml-dev` may have to be installed (required for `pbio`):
   
    
```
# Optional: update gcc if necessary.
conda install -c salford_systems libgcc-6=6.2.0

# Optional: install llvm
conda install -c anaconda llvm
conda install -c numba llvmlite

# Create the anaconda environment.
conda create -n my_new_environment python=3.6 anaconda

# Activate the new environment.
source activate my_new_environment

# Clone the git repository
git clone git@github.com:dieterich-lab/rp-bp.git
cd rp-bp

pip --verbose install -r requirements.txt
```

[Back to top](#toc)

<a id='uninstallation'></a>

## Uninstallation

The `rpbp` package and related dependencies can be removed with pip:

``pip uninstall pbio rpbp``

If the packages were installed in a dedicated virtual environment, this environment can also be cleared or removed.