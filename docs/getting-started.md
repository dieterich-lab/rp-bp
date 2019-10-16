# Getting started

## Ribosome profiling with Bayesian predictions

Source code available at [https://github.com/dieterich-lab/rp-bp](https://github.com/dieterich-lab/rp-bp)

### Installation

This package is written in Python3. It has a number of external dependencies, mostly standard bioinformatics tools. Rp-Bp is not published on PyPI, but the installation is easily managed through `pip3`. The required privileges are determined by the installation location of `pip3`. In particular, if `pip3` does not require sudo access, then none of the installation process requires sudo access; this is the case within a virtual environment, for example. For detailed instructions, including dependencies/prerequisites and step-by-step details of installing within a virtual environment and anaconda. refer to [installation](installation.html).

### Usage

Please see [Running the Rp-Bp pipeline step-by-step](usage-instructions.html) for more detailed usage instructions. We also provide a number of tools to "post-process" and visualise the results, see [QC and downstream analysis of the Rp-Bp results](analysis-scripts.html). To get started, the package also includes a small example using a *C. elegans* dataset. Please see [Running Rp-Bp on the example dataset](running-example.html) for instructions on running the example.

### Citation

Malone, B.; Atanassov, I.; Aeschimann, F.; Li, X. & Dieterich, C. Bayesian prediction of RNA translation from ribosome profiling. Nucleic Acids Research, 2017, gkw1350. (Volume and pages have not yet been assigned). The paper is [available at NAR](https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkw1350).

### License

The MIT License (MIT)

Copyright (c) 2016 dieterich-lab

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## How to report issues

To report bugs and issues with the `rpbp` package, please use the [bug tracker](https://github.com/dieterich-lab/rp-bp/issues). Please
follow the instructions and guidelines given in the template.
