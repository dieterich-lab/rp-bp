
# Installing the Rp-Bp pipeline

This document describes how to install the `rpbp` package. 

The `rpbp` package was developed and tested on Debian-based systems, but should also run with any distribution.

<a id="toc"></a>

* [Prerequisites](#prerequisites)
* [Installation](#installation)
    * [Virtual environment installation](#virtual-environment-installation)
    * [Anaconda installation](#anaconda-installation)
    * [Using docker and a precompiled container](#Using-docker-and-a-precompiled-container) 
* [Uninstallation](#uninstallation)

---

<a id="prerequisites"></a>

## Prerequisites

The `rpbp` package requires the following external dependencies:

* [Flexbar](https://github.com/seqan/flexbar), version 3.5.0
* [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), version 2.3.0
* [STAR](https://github.com/alexdobin/STAR), version 2.6.1d
* [Samtools](http://www.htslib.org), version 1.7

These dependencies are *NOT* installed automatically, and must be installed prior to installing
the `rpbp` package. Do not forget to add these programs to your `$PATH` or to symlink the executable/binary files
to your `~/.local/bin` directory! Alternatively, if your are installing `rpbp` in 
an [anaconda virtual environment](#anaconda-installation), these dependencies (except for STAR) can
be installed with `conda`.

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

To install the local VCS project in development mode, use the `--editable` or `-e` option, otherwise
this flag can be ignored. 

The `--force-recompile` flag can also be passed via `install-option` to force the recompilation of the Stan models.
This flag does *NOT* normally have to be set, unless one has to 
re-install the `rpbp` package (*e.g* after upgrade or changing versions of Pystan). Note that 
the use of `--install-option` currently disables all use of wheels.

** Note: `install-option` seems to leak across lines... Until further testing is done,
use of the `--force-recompile` flag is not recommended.

Pinned version of selected dependencies are installed with the `pbio` package via
the `requirements.txt` file for reproducible installation.

<a id='virtual-environment-installation'></a>

## Installation (virtual environment)

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
pip --verbose install -r requirements.txt [-e] . [--install-option="--force-recompile"] 2>&1 | tee install.log

```

[Back to top](#toc)

<a id='anaconda-installation'></a>

## Anaconda installation

The package can also be installed within an [anaconda](https://www.continuum.io/) environment. 

   * The version of gcc included by default may not properly compile the Stan model files, so it must be updated first.
   * In addition, `llvm` may have to be installed (required for `pbio`):
    
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
pip --verbose install -r requirements.txt [-e] . [--install-option="--force-recompile"] 2>&1 | tee install.log

```

[Back to top](#toc)

<a id='uninstallation'></a>

## Using docker and a precompiled container

To run the `rpbp` package and its dependencies in an isolated container, you first need to install docker. For this please follow the instructions on the [docker](https://docker.com) webpage. 

To interactively run and create an executable [docker](https://docker.com) container layer, you have to run:
```
docker run -it --rm --mount type=bind,source=/path/to/rpbp/input/and/reference,target=/target/mount/path/in/docker/container dquz/rp-bp
```
Further descriptions on how to create and execute docker containers can be found at [docker](https://docker.com).

[Back to top](#toc)

<a id='anaconda-installation'></a>

## Uninstallation

The `rpbp` package and related Python dependencies can be removed with pip:

``pip uninstall pbio rpbp``

If the packages were installed in a dedicated virtual environment, this environment can also be cleared or removed.
