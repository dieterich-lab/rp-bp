This project contains scripts for running the Rp-Bp and Rp-chi translation prediction pipelines.

#Installation instructions

This package is written in python3. Internally, pip3 is used to install the python dependencies. It has a number of external dependencies, mostly standard bioinformatics tools. Please see [docs/usage-instructions.md](docs/usage-instructions.md#prerequisites). (Locally, ``docs/usage-instructions.html`` may be easier to view.)

Installation is managed through pip3. The required privileges are determined by the installation location of pip3. In particular, if pip3 does not require sudo access, then none of the installation process requires sudo access; this is the case within a virtual environment, for example.

Installation requires running make:

``pip3 install --verbose .``


If possible, we recommend installing inside a virtual environment.

Please see [docs/installation.md](docs/installation.md) for more detailed installation instructions, including prerequisites and step-by-step details of installing within a virtual machine. (Locally, ``docs/installation.html`` may be easier to view.)

#Uninstallation instruction

The rp-bp python packages can also be removed with pip.

``pip3 uninstall misc riboutils rpbp``


#Usage instructions

Please see [docs/usage-instructions.md](docs/usage-instructions.md) for usage. (Locally, ``docs/usage-instructions.html`` may be easier to view.)

#Example instructions

The package includes a small example using a *C. elegans* dataset. Please see [docs/running-example.md](docs/running-example.md) for instructions on running the example. (Locally, ``docs/running-example.html`` may be easier to view.)
