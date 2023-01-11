Getting started
###############

**Rp-Bp** can be used for ORF discovery, or simply to estimate periodicity in a set of Ribo-seq replicates, to know which samples and read lengths are usable for downstream analyses.

To get started, you need *(i)* Ribo-seq data (FASTQ), *(ii)* genome sequence and annotation for your organism (FASTA, GTF), and *(iii)* ribosomal sequences for *in-silico* rRNA removal (FASTA). We also recommend to filter protocol-specific or general adapter sequences, this can be done by providing a FASTA file containing a set of adapter sequences.

.. _getting_started_containers:

Containers
**********

The quickly get started, the preferred way is to use a Docker container with **Rp-Bp** pre-installed:

.. code-block:: bash
    docker pull ...


.. _getting_started_install:

Installation
************

Supported operating systems
---------------------------

Debian GNU/Linux and macOS.


To install from Bioconda using Conda, set up `Bioconda <https://bioconda.github.io/#usage>`_ or `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ if not already done, then:

.. code-block:: bash
    conda install rpbp

The Conda installation will resolve all package dependencies.


**Rp-Bp** is also available on `PyPI <https://pypi.org/project/rpbp>`_, but required dependencies are not installed automatically, and must be installed prior to running the pipeline.


For detailed instructions, refer to `Installation <installation.html>`_.


General usage
*************

For detailed usage instructions, refer to `Running the pipeline <usage.html>`_.
For visualisation and quality control, see `Visualisation and quality control <apps.html>`_.

To get started, the package also includes a small example using a *C. elegans* dataset.
To run the example, refer to `Tutorial <tutorial.html>`_.


How to report issues
====================

Bugs and issues should be reported in the `bug tracker <https://github.com/dieterich-lab/rp-bp/issues>`_. Follow the instructions and guidelines given in the template.


How to cite
===========

Brandon Malone, Ilian Atanassov, Florian Aeschimann, Xinping Li, Helge Großhans, Christoph Dieterich. `Bayesian prediction of RNA translation from ribosome profiling <https://doi.org/10.1093/nar/gkw1350>`_, *Nucleic Acids Research*, Volume 45, Issue 6, 7 April 2017, Pages 2960-2972.


License
=======

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
