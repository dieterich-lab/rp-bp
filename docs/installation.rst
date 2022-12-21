Installation
############

In general, we recommend to use `containers <getting-started.html#Containers>`_, however if you want to install **Rp-Bp**, this can be done with Conda (preferred) or PyPI. 


.. _conda_install:

Conda installation
******************



.. _pypi_install:

PyPI installation
*****************

If you already have the required dependencies installed on your system, to install from `PyPI <https://pypi.org/project/rpbp>`_

.. code-block:: bash
    # create a virtual environment...
    python3 -m venv /path/to/virtual/environment 
    # ... activate it ... 
    source /path/to/virtual/environment/bin/activate
    # ... and install Rp-Bp
    pip install rpbp
    
To install the local VCS project in development mode, use the ``--editable`` or ``-e`` option, otherwise this flag can be ignored.

Required dependencies
---------------------

* `Flexbar <https://github.com/seqan/flexbar>`_
* `Bowtie 2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_
* `STAR <https://github.com/alexdobin/STAR>`_
* `Samtools <http://www.htslib.org>`_

.. important:: 
    Required dependencies are NOT installed automatically, and must be installed prior to running the pipeline. You may have to add these programs to your ``$PATH`` or install pre-compiled binaries under *e.g.* ``/home/$USER/.local/bin``.


.. _uninstall:
    
Uninstallation
**************

To remove **Rp-Bp** if installed with pip:

.. code-block:: bash
    pip uninstall rpbp

To remove **Rp-Bp** if installed with conda:

.. code-block:: bash
    conda uninstall rpbp

If the package is installed in a dedicated virtual (conda) environment, this environment can also be cleared or removed.
