# **R**ibosome **p**rofiling with **B**ayesian **p**redictions (Rp-Bp)

An unsupervised Bayesian approach to predict translated open reading frames (ORFs) from ribosome profiles, using an automatic **B**ayesian **P**eriodic fragment length and ribosome **P**-site offset **S**election (BPPS).

![rpbp](docs/images/logo-rpbp-final.png)

---

## Documentation

Read the Docs: [![Documentation Status](https://readthedocs.org/projects/rp-bp/badge/?version=latest)](http://rp-bp.readthedocs.io/en/latest/?badge=latest) or just click [here](http://rp-bp.readthedocs.io/en/latest/) to access the complete documentation.

## Installation

This package is written in Python3. It has a number of external dependencies, mostly standard bioinformatics tools. Rp-Bp is not published on PyPI, but the installation is easily managed through `pip3`. The required privileges are determined by the installation location of `pip3`. In particular, if `pip3` does not require sudo access, then none of the installation process requires sudo access; this is the case within a virtual environment, for example. For detailed instructions, including dependencies/prerequisites and step-by-step details of installing within a virtual environment and anaconda. refer to [installation instructrions](http://rp-bp.readthedocs.io/en/latest/installation.html).

<a name="get-start-usage"></a>

## Usage

Please see [Running the Rp-Bp pipeline step-by-step](http://rp-bp.readthedocs.io/en/latest/usage-instructions.html) for more detailed usage instructions. We also provide a number of tools to "post-process" and visualise the results, see [QC and downstream analysis of the Rp-Bp results](http://rp-bp.readthedocs.io/en/latest/analysis-scripts.html). To get started, the package also includes a small example using a _C. elegans_ dataset. Please see [Running Rp-Bp on the example dataset](http://rp-bp.readthedocs.io/en/latest/running-example.html) for instructions on running the example.

<a name="get-start-cite"></a>

## Citation

Malone, B.; Atanassov, I.; Aeschimann, F.; Li, X. & Dieterich, C. Bayesian prediction of RNA translation from ribosome profiling. Nucleic Acids Research, 2017, gkw1350. (Volume and pages have not yet been assigned). The paper is [available at NAR](https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkw1350).
