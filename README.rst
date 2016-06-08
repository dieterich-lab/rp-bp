This project contains scripts for running the Rp-Bp and Rp-chi translation prediction pipelines.

**Install**

This package is written in python3. Internally, pip3 is used to install the python dependencies. Additionally, the `OpenBLAS <http://www.openblas.net/>`_ library is used for efficiency. The installation procedure also downloads and compiles this library. Automatic installation of OpenBLAS has been thoroughly tested within a virtual environment (see below) on Ubuntu. Other environments may require adjustments. Please contact us for help.

Installation is managed through the included makefile. The required privileges are determined by the installation location of pip3. In particular, if pip3 does not require sudo access, then none of the installation process requires sudo access; this is the case within a virtual environment, for example.

Installation simply requires running make:

``make``

The python packages and OpenBlas installation can also be removed with make:

``make clean``

If possible, we recommend installing inside a virtual environment. See `here 
<http://www.simononsoftware.com/virtualenv-tutorial-part-2/>`_, for example. The virtualenvwrapper package seems to work better in Python2, but python3 must be used inside the virtual environment.

**Complete installation instructions within a virtual environment**

These instructions explain how to install the software and most dependencies from scratch without required root access.
It only requires standard development libraries and tools, like gcc and the gzip development headers.
The python build scripts will also output a line like "The necessary bits to build these optional modules were not found" if any optional libraries, developoment headers, etc., are not found.

The commands below are presumably executed in a directory like ``$HOME/install``.
They install the python executables into ``$HOME/local/bin`` (or wherever the ``prefix`` option is located).
So that directory must be in the ``$PATH``.
This can be accomplished by adding a line like ``export PATH=$HOME/local/bin:$PATH`` in the file ``.bashrc`` on Ubuntu.

* Download, extract and install Python 2. ``wget https://www.python.org/ftp/python/2.7.11/Python-2.7.11.tgz && tar -xvf Python-2.7.11.tgz && cd Python-2.7.11 && ./configure --prefix=$HOME/local --with-ensurepip=upgrade && make && make install && cd ..``
* Download, extract and install Python 3. ``wget https://www.python.org/ftp/python/3.5.1/Python-3.5.1.tgz && tar -xvf Python-3.5.1.tgz && cd Python-3.5.1 && ./configure --prefix=$HOME/local --with-ensurepip=upgrade && make && make install && cd ..``
* Upgrade both versions of pip. ``pip2 install --upgrade pip && pip3 install --upgrade pip``
* Install the virtual environment wrapper for Python 2. ``pip2 install --upgrade virtualenvwrapper``
* Add the required lines to .bashrc: ``export WORKON_HOME=$HOME/.virtualenvs`` and ``source $HOME/local/bin/virtualenvwrapper_lazy.sh``. See `here <http://www.simononsoftware.com/virtualenv-tutorial-part-2/>`_ for more details.
* Add the lines required for the OpenBLAS library to .bashrc: ``export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH``
* Create a virtual environment using Python 3. ``mkvirtualenv rpbp -p $HOME/local/bin/python3``
* Clone the git repository. ``git clone git@github.com:dieterich-lab/rp-bp.git``
* Change into the ``rp-bp`` directory and build the package. ``cd rp-bp && make``

The build process includes compiling several libraries for optimized numerical calculations. Due to the optimized nature of these libraries, the initial installation can take up to an hour.

To use the programs in the future, ensure the virtual environment is active. ``workon rpbp``.

**Installation instructions without using a virtual environment**

We recommend installing the application in a virtual environment as described above.
If this is not desired for some reason, the following instructions can be used to install the package without sudo access in a user's home directory.

The commands below are presumably executed in a directory like ``$HOME/install``.
They install the python executables into ``$HOME/local/bin`` (or wherever the ``prefix`` option is located).
So that directory must be in the ``$PATH``.
This can be accomplished by adding a line like ``export PATH=$HOME/local/bin:$PATH`` in the file ``.bashrc`` on Ubuntu.

* Download, extract and install Python 3. ``wget https://www.python.org/ftp/python/3.5.1/Python-3.5.1.tgz && tar -xvf Python-3.5.1.tgz && cd Python-3.5.1 && ./configure --prefix=$HOME/local --with-ensurepip=upgrade && make && make install && cd ..``
* Upgrade pip. ``pip3 install --upgrade pip``
* Clone the git repository. ``git clone git@github.com:dieterich-lab/rp-bp.git``
* Change into the ``rp-bp`` directory and build the package. ``cd rp-bp && make``

The build process includes compiling several libraries for optimized numerical calculations. Due to the optimized nature of these libraries, the initial installation can take up to an hour.

**Usage instructions**

Please see the usage-instructions.html file for usage.
