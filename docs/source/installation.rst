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

There is no *latest* tag, you need to specify the version tag. See `rpbp/tags <https://quay.io/repository/biocontainers/rpbp?tab=tags>`_ for valid values for <tag>.


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

.. _pypi_install:

Contributing to **Rp-Bp**
-------------------------

To install the local VCS project in development mode

.. code-block:: bash

    # create a conda environment...
    conda create -n rpbp_dev
    # ...activate it...
    conda activate rpbp_dev
    # ... and only install dependencies
    (rpbp_dev) conda install --only-deps rpbp
    # clone the git repository
    (rpbp_dev) git clone https://github.com/dieterich-lab/rp-bp.git && cd rp-bp
    # install
    (rpbp_dev) pip --verbose install --no-deps -e . 2>&1 | tee install.log


PyPI installation
^^^^^^^^^^^^^^^^^

We do not recommend to install **Rp-Bp** directly from `PyPI <https://pypi.org/project/rpbp>`_.
However, if you already have the required dependencies installed on your system, to install

.. code-block:: bash

    # create a virtual environment...
    python3 -m venv rpbp_pypi
    # ... activate it ...
    source rpbp_pypi/bin/activate
    # ... and install Rp-Bp
    (rpbp_pypi) pip install rpbp


Required dependencies
"""""""""""""""""""""

* `Flexbar <https://github.com/seqan/flexbar>`_
* `Bowtie 2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_
* `STAR <https://github.com/alexdobin/STAR>`_
* `Samtools <http://www.htslib.org>`_

.. warning::

    Conda installation or containers include all dependencies. With a PyPI installation, you need to install required dependencies. Executables or binaries must be in your ``$PATH``.

.. _uninstall:

Uninstallation
--------------

Remove the conda environment

.. code-block:: bash

    conda env remove --name rpbp

or remove the package installed in another environment

.. code-block:: bash

    # remove the rpbp package from myenv environment...
    (myenv) conda remove -n myenv rpbp


To remove **Rp-Bp** if installed with pip

.. code-block:: bash

    pip uninstall rpbp


If the package is installed in a dedicated python virtual environment, this environment can also be removed.
