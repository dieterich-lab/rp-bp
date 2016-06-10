This project contains scripts for running the Rp-Bp and Rp-chi translation prediction pipelines.

**Install**

This package is written in python3. Internally, pip3 is used to install the python dependencies. Additionally, the `OpenBLAS <http://www.openblas.net/>`_ library is used for efficiency. The installation procedure also downloads and compiles this library. Automatic installation of OpenBLAS has been thoroughly tested within a virtual environment (see below) on Ubuntu. Other environments may require adjustments. Please contact us for help.

Installation is managed through the included makefile. The required privileges are determined by the installation location of pip3. In particular, if pip3 does not require sudo access, then none of the installation process requires sudo access; this is the case within a virtual environment, for example.

Installation simply requires running make:

``make``

The python packages and OpenBlas installation can also be removed with make:

``make clean``

If possible, we recommend installing inside a virtual environment.

Please see [installation.html](installation.html) file for more detailed installation instructions, including prerequisites and step-by-step details of installing within a virtual machine.

**Usage instructions**

Please see the [usage-instructions.html](usage-instructions.html) file for usage.
