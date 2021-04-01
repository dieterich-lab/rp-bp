# Getting started

## Ribosome profiling with Bayesian predictions

We propose an unsupervised Bayesian approach to predict translated open reading frames (ORFs) from ribosome profiles. We rely on MCMC sampling using `Stan` to estimate posterior distributions
of the likelihood of translation for each identified ORF. We implement an automatic Bayesian Periodic fragment length and ribosome P-site offset Selection (BPPS), *i.e.* read lengths and ribosome P-site offsets are inferred from the data, without supervision. Hence our method is able to handle *de novo* translatome annotation by directly assessing the periodicity of the Ribo-seq signal.

### Supported operating systems

Our workflow was developed and tested on Debian GNU/Linux, including the Ubuntu distribution.
macOS is currently not fully supported.

### Installation

`rpbp` is written in Python3. The package is not currently published on PyPI, but the installation
is done through `pip`, after cloning the Git repository. Note that `pip` will not install external dependencies
that are required to run the pipeline. Supported versions of Python and major dependencies are specified in each release.

For detailed instructions, refer to [installation](installation.html).

### Usage

For detailed usage instructions, refer to [Running the Rp-Bp pipeline step-by-step](usage-instructions.html).
For quality control and visualisation of the results, see [QC and downstream analysis of the Rp-Bp results](analysis-scripts.html).

To get started, the package also includes a small example using a *C. elegans* dataset.
To run the example, refer to [Running Rp-Bp on the example dataset](running-example.html).

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

## How to cite

Brandon Malone, Ilian Atanassov, Florian Aeschimann, Xinping Li, Helge Großhans, Christoph Dieterich.
[Bayesian prediction of RNA translation from ribosome profiling](https://doi.org/10.1093/nar/gkw1350), *Nucleic Acids Research*, Volume 45, Issue 6, 7 April 2017, Pages 2960-2972.


