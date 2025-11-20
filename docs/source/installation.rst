.. _installation_full:

Installation
============

Containers
----------

To use a container (Docker or Singularity) with **Rp-Bp** pre-installed, simply pull, and you're done!

.. code-block:: bash

    # docker or...
    docker pull quay.io/biocontainers/rpbp:<tag>
    # ...singularity
    singularity pull rpbp.sif docker://quay.io/biocontainers/rpbp:<tag>

There is no *latest* tag, you need to specify the version tag. See `rpbp/tags <https://quay.io/repository/biocontainers/rpbp?tab=tags>`_ for valid values for <tag>. Check the `Tutorials <tutorial.html>`_ on how to use the containers.


.. _conda_install:

Conda installation
------------------

If required, set up the conda channels as described `here <https://bioconda.github.io/#usage>`_, and install with

.. code-block:: bash

    # preferably install in some conda environment...
    conda install rpbp

or create an environment, called *rpbp*, containing the **Rp-Bp** package

.. code-block:: bash

    conda create -n rpbp rpbp

.. tip::

    `Mamba <https://mamba.readthedocs.io/en/latest/user_guide/mamba.html#mamba>`_ can be used as a drop-in replacement, you can swap almost all commands between conda and mamba.

.. _pypi_install:

Contributing to **Rp-Bp**
-------------------------

To install the local VCS project in development mode

.. code-block:: bash

    # create a new environment...
    mamba create -n rpbp
    # ... activate...
    mamba activate rpbp
    # ... and install dependencies...
    mamba install --only-deps rpbp
    # ... clone the git repository and install the package
    git clone https://github.com/dieterich-lab/rp-bp.git && cd rp-bp
    pip install --no-deps --editable . 2>&1 | tee install.log

Alternatively, clone the git repository and install from the `yaml spec file <https://mamba.readthedocs.io/en/latest/user_guide/mamba.html#conda-yaml-spec-files>`_

.. code-block:: bash

   mamba env create -f environment.yml
   # ... activate environment...
   mamba activate rpbp
   pip install --no-deps --editable . 2>&1 | tee install.log

Finally, install test dependencies.

PyPI installation
^^^^^^^^^^^^^^^^^

To install the package from `PyPI <https://pypi.org/project/rpbp>`_

.. code-block:: bash

    # create a virtual environment...
    python3 -m venv rpbp
    # ... activate ...
    source rpbp/bin/activate
    # ... and install the package
    pip install rpbp

**Required dependencies:** `Flexbar <https://github.com/seqan/flexbar>`_, `Bowtie 2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_, `STAR <https://github.com/alexdobin/STAR>`_, `Samtools <http://www.htslib.org>`_, and `FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_.

.. warning::

    A PyPI installation only installs the python package. You need to install required dependencies separately. Executables or binaries must be in your ``$PATH``.

.. _uninstall:

Uninstallation
--------------

Remove the environment

.. code-block:: bash

    mamba env remove --name rpbp

or remove the package installed in another environment

.. code-block:: bash

    # remove the rpbp package from myenv environment...
    mamba remove -n myenv rpbp

To remove **Rp-Bp** if installed with pip

.. code-block:: bash

    pip uninstall rpbp

If the package is installed in a dedicated python virtual environment, remove this environment.
