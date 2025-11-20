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
    # create a conda environment called rpbp and install rpbp
    conda create -n rpbp rpbp

or use a container

.. code-block:: bash

    # docker or...
    docker pull quay.io/biocontainers/rpbp:<tag>
    # ...singularity
    singularity pull rpbp.sif docker://quay.io/biocontainers/rpbp:<tag>

There is no *latest* tag, you need to specify the version tag. See `rpbp/tags <https://quay.io/repository/biocontainers/rpbp?tab=tags>`_ for valid values for ``<tag>``.

For detailed installation instructions, refer to :ref:`installation_full`.


**Rp-Bp** quickstart
--------------------

In a nutshell, you need to prepare genome indices and annotations for your organism by calling

.. code-block:: bash

    prepare-rpbp-genome [options] config


To estimate periodicity on a set of Ribo-seq samples, or to run the ORF discovery pipeline, simply call

.. code-block:: bash

    run-all-rpbp-instances [options] config

For more information and guidelines how to prepare the configuration file and run the pipeline, refer to the :ref:`user_guide`.
For visualization and quality control, see :ref:`apps`.

To get started, the package includes a small example dataset. Check the :ref:`all_tutorials`.

How to report issues
--------------------

Bugs and issues should be reported in the `bug tracker <https://github.com/dieterich-lab/rp-bp/issues>`_. Follow the instructions and guidelines given in the template.


How to contribute
-----------------

Contributions are welcome! New code should follow `Black <https://black.readthedocs.io/en/stable/>`_ and `flake8 <https://flake8.pycqa.org/en/latest/>`_. Install development dependencies inside a virtual environment, see :ref:`pypi_install`. A typical development workflow would include *(i)* forking the repository, *(ii)* creating a new branch for your PR, *(iii)* adding features or bug fixes, *(iv)* making sure all tests are passing, *(v)* building the documentation if necessary, and *(vi)* opening a PR back to the main repository. If you're fixing a bug, add a test. Run it first to confirm it fails, then fix the bug, and run it again to confirm it's fixed. If adding a new feature, add a test, or first open an issue to discuss the idea.

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
