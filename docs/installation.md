
# Installing the Rp-Bp pipeline

This document describes how to install the `rpbp` package. 

<a id="toc"></a>

* [Supported operating systems](#os)
* [Required dependencies](#prerequisites)
* [Installation](#installation)
    * [Virtual environment installation](#virtual-environment-installation)
    * [Anaconda installation](#anaconda-installation)
* [Uninstallation](#uninstallation)

---

<a id="os"></a>

## Supported operating systems

Our workflow was developed and tested on Debian GNU/Linux (Debian jessie 8 64-bit), including the Ubuntu distribution (Ubuntu Bionic Beaver 18.04.3 LTS 64-bit).
macOS is currently not fully supported.

<a id="prerequisites"></a>

## Required dependencies

### External software

* [Flexbar [== 3.5.0]](https://github.com/seqan/flexbar)
* [Bowtie 2 [== 2.3.0]](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [STAR [== 2.6.1d]](https://github.com/alexdobin/STAR)
* [Samtools [== 1.7]](http://www.htslib.org)

!!! important
    These dependencies are NOT installed automatically, and must be installed prior to running the `rpbp` pipeline. You may have to add these programs to your `$PATH` or install pre-compiled binaries under `/home/$USER/.local/bin`. Alternatively, if your are installing `rpbp` in an [anaconda virtual environment](#anaconda-installation), these dependencies (except for STAR) can be installed with `conda`.

### Python packages

Supported versions of Python and major Python dependencies are specified in each
release, and installed automatically. To avoid conflicting versions, we recommend
installing `rpbp` in a virtual environment.

- `rpbp` version 2.0.0 is tested for Python \[>=3.6,<3.7.0a0\]

!!! note
    If installation fails due to missing `OpenBLAS` dependencies for `scipy`, please follow the instructions [here](https://gist.github.com/bmmalone/1b5f9ff72754c7d4b313c0b044c42684).

[Back to top](#toc)

<a id='installation'></a>

## Installation

To install the local VCS project in development mode, use the `--editable` or `-e` option, otherwise
this flag can be ignored. The GitHub installation will install the most recent version directly from the source repository.

The `--force-recompile` flag can also be passed via `install-option` to force the recompilation of the Stan models.
This flag does *NOT* normally have to be set, unless one has to 
re-install the `rpbp` package (*e.g* after upgrade or changing versions of Pystan). Note that 
the use of `--install-option` currently disables all use of wheels.

!!! warning
    `install-option` seems to leak across lines... Until further testing is done, use of the `--force-recompile` flag is not recommended!

Pinned version of selected dependencies are installed with the `pbio` package via
the `requirements.txt` file for reproducible installation.

<a id='virtual-environment-installation'></a>

### Installation (virtual environment, recommended)

To install the `rpbp` package and its dependencies, first create a virtual environment:
 
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
git clone https://github.com/dieterich-lab/rp-bp.git
cd rp-bp

# The period is required, it is the local project path (rp-bp)
pip --verbose install -r requirements.txt . [--install-option="--force-recompile"] 2>&1 | tee install.log

```

[Back to top](#toc)

<a id='anaconda-installation'></a>

### Anaconda installation

The package can also be installed within a [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) virtual environment.

!!! note
    The version of `gcc` included by default may not properly compile the Stan model files, so it must be updated first. In addition, `llvm` may have to be installed (required for `pbio`).

```
# Optional: update gcc and install llvm.
conda install -c salford_systems libgcc-6=6.2.0
conda install -c anaconda llvm
conda install -c numba llvmlite

# Create the anaconda environment.
conda create -n my_new_environment python=3.6 anaconda

# Activate the new environment.
source activate my_new_environment

# Optional: install external dependencies (except STAR)
conda install -c bioconda flexbar=3.5.0
conda install -c bioconda bowtie2=2.3.0
conda install -c bioconda samtools=1.7

# Clone the git repository
git clone https://github.com/dieterich-lab/rp-bp.git
cd rp-bp

# The period is required, it is the local project path (rp-bp)
pip --verbose install -r requirements.txt . [--install-option="--force-recompile"] 2>&1 | tee install.log

```

[Back to top](#toc)

<a id='uninstallation'></a>

## Uninstallation

The `rpbp` package and related Python dependencies can be removed with pip:

``pip uninstall pbio rpbp``

If the packages were installed in a dedicated virtual environment, this environment can also be cleared or removed.
