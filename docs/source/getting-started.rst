Getting started
===============

What is **Rp-Bp**?
------------------

**Rp-Bp** is an unsupervised Bayesian approach to predict translated open reading frames (ORFs) from ribosome profiles. **Rp-Bp** can be used for ORF discovery, or simply to estimate periodicity in a set of Ribo-seq samples.

To get started, you need

* Ribo-seq data (FASTQ)
* genome sequence and annotation for your organism (FASTA, GTF)
* ribosomal sequence for *in-silico* rRNA removal (FASTA)
* protocol-specific or general adapter sequences to be removed (FASTA)

.. _getting_started:

Installation
------------

Install with

.. code-block:: bash

    # set up the conda channels if required
    # conda config --add channels defaults
    # conda config --add channels bioconda
    # conda config --add channels conda-forge
    # conda config --set channel_priority strict

    # create a conda environment called rpbp and install rpbp
    conda create -n rpbp rpbp

or use a container

.. code-block:: bash

    # docker or...
    docker pull quay.io/biocontainers/rpbp:<tag>
    # ...singularity
    singularity pull rpbp.sif docker://quay.io/biocontainers/rpbp:<tag>

There is no *latest* tag, you need to specify the version tag. See `rpbp/tags <https://quay.io/repository/biocontainers/rpbp?tab=tags>`_ for valid values for ``<tag>``.

For detailed installation instructions, refer to `Installation <installation.html>`_.


**Rp-Bp** quickstart
--------------------

In a nutshell, you need to prepare genome indices and annotations for your organism by calling

.. code-block:: bash

    prepare-rpbp-genome <config> [options]


To estimate periodicity on a set of Ribo-seq samples, or to run the ORF discovery pipeline, simply call

.. code-block:: bash

    run-all-rpbp-instances <config> [options]

To get started, the package also includes a small example using a *C. elegans* dataset. Check the `Tutorials <tutorial.html>`_.

For more information and guidelines on how to prepare the configuration file, refer to the `User guide <user-guide.html>`_.
For visualisation and quality control, see `Visualization and quality control <apps.html>`_.


How to report issues
--------------------

Bugs and issues should be reported in the `bug tracker <https://github.com/dieterich-lab/rp-bp/issues>`_. Follow the instructions and guidelines given in the template.


How to contribute
-----------------

Contributions are welcome! New code should follow `Black <https://black.readthedocs.io/en/stable/>`_ and `flake8 <https://flake8.pycqa.org/en/latest/>`_. A typical development workflow would include *(i)* forking the repository, *(ii)* creating a new branch for your PR, *(iii)* adding features or bug fixes, *(iv)* making sure all tests are passing, *(v)* building the documentation if necessary, and *(vi)* opening a PR back to the main repository.

Running the tests
^^^^^^^^^^^^^^^^^

We use pytest to test **Rp-Bp**. Currently, only regression tests are implemented. Dependencies can be installed with ``pip install -e .[tests]``.

Building the docs
^^^^^^^^^^^^^^^^^

Dependencies for building the documentation can be installed with ``pip install -e .[docs]``.

Semantic versioning
^^^^^^^^^^^^^^^^^^^

We try to follow `semantic versioning <https://semver.org/>`_.


How to cite
-----------

Brandon Malone, Ilian Atanassov, Florian Aeschimann, Xinping Li, Helge Großhans, Christoph Dieterich. `Bayesian prediction of RNA translation from ribosome profiling <https://doi.org/10.1093/nar/gkw1350>`_, *Nucleic Acids Research*, Volume 45, Issue 6, 7 April 2017, Pages 2960-2972.


License
-------

The MIT License (MIT). Copyright (c) 2016 dieterich-lab.
