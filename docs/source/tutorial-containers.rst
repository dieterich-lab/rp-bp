.. _tutorial_containers:

How to use the containers
=========================

First pull a Docker or Singularity container, see `installation <installation.html>`_. For this tutorial, we use a general mechanism for persisting data, which allows to create and modify files on the host system from within the container. We use the *c-elegans-chrI-example* with the *c-elegans-test.yaml* configuration file, see also `How to run the pipeline <tutorial-cel.html>`_.

.. note::

    You can also launch a container with an interactive shell *e.g.* ``docker run -it quay.io/biocontainers/rpbp:3.0.1--py310h30d9df9_0 bash`` or ``singularity shell rpbp.sif``. With ``singularity shell``, ``$HOME`` is mounted by default.

.. attention::

    In the following, do not forget to modify the tag ``3.0.1--py310h30d9df9_0`` according to what you pulled! For Singularity, adjust the name of the Singularity image format file ``rpbp.sif`` and/or the path according to your needs.


How to run the pipeline
-----------------------

To run the pipeline, change the paths in the configuration file to point to the location where the directory is mounted in the container. You can do this using a text editor, or simply by modifying the file in place

.. code-block:: bash

    sed -i 's|/path/to/your/c-elegans-example|/data|g' c-elegans-test.yaml


.. important::

    Default parameters were modified for the example and included in the configuration file. If you use this configuration file as a general template for your data, do not forget to remove everything below the line "REMOVE BELOW THIS LINE IF YOU USE THIS CONFIGURATION FILE AS TEMPLATE FOR YOUR DATA".


You can now create the genome indices and annotations. For this small example, it is important to downscale the STAR ``--genomeSAindexNbases`` parameter as follows


.. code-block:: bash

    docker run --volume `pwd`:/data quay.io/biocontainers/rpbp:3.0.1--py310h30d9df9_0 prepare-rpbp-genome /data/c-elegans-test.yaml --star-options "--genomeSAindexNbases 10" --num-cpus 4 --logging-level INFO --log-file /data/rpbp-genome.log


.. code-block:: bash

    singularity run --bind `pwd`:/data rpbp.sif prepare-rpbp-genome /data/c-elegans-test.yaml --star-options "--genomeSAindexNbases 10" --num-cpus 4 --logging-level INFO --log-file /data/rpbp-genome.log


The file *rpbp-genome.log* contains logging output for the reference preprocessing. You now have a new directory called *WBcel235.79.chrI* with genome indices and annotations.

Finally, run the ORF discovery pipeline


.. code-block:: bash

    docker run --volume `pwd`:/data quay.io/biocontainers/rpbp:3.0.1--py310h30d9df9_0 run-all-rpbp-instances /data/c-elegans-test.yaml --merge-replicates --run-replicates --keep-intermediate-files --num-cpus 4 --logging-level INFO --log-file /data/rpbp-pipeline.log


.. code-block:: bash

    singularity run --bind `pwd`:/data rpbp.sif run-all-rpbp-instances /data/c-elegans-test.yaml --merge-replicates --run-replicates --keep-intermediate-files --num-cpus 4 --logging-level INFO --log-file /data/rpbp-pipeline.log


The file *rpbp-pipeline.log* contains logging output for the different processing steps. You now have four new directories (*with-*, *without-*) including output from Flexbar, Bowtie2, and STAR, and directories with **Rp-Bp** output: *metagene-profiles*, *orf-profiles*, and *orf-predictions*. The *orf-predictions* include the output for each sample *c-elegans-rep-1* and *c-elegans-rep-2* as well as for the merged replicates *c-elegans-test*.

How to summarize the results and launch the apps
------------------------------------------------

Prepare the summary output for the *profile construction dashboard*

.. code-block:: bash

    docker run --volume `pwd`:/data quay.io/biocontainers/rpbp:3.0.1--py310h30d9df9_0 summarize-rpbp-profile-construction /data/c-elegans-test.yaml --num-cpus 4 --logging-level INFO --log-file /data/rpbp-profile-summary.log


.. code-block:: bash

    singularity run --bind `pwd`:/data rpbp.sif summarize-rpbp-profile-construction /data/c-elegans-test.yaml --num-cpus 4 --logging-level INFO --log-file /data/rpbp-profile-summary.log


and for the *predictions dashboard*

.. code-block:: bash

    docker run --volume `pwd`:/data quay.io/biocontainers/rpbp:3.0.1--py310h30d9df9_0 summarize-rpbp-predictions /data/c-elegans-test.yaml --no-replicates --circos-bin-width 10000 --circos-show-chroms I --logging-level INFO --log-file /data/rpbp-predictions-summary.log

.. code-block:: bash

    singularity run --bind `pwd`:/data rpbp.sif summarize-rpbp-predictions /data/c-elegans-test.yaml --no-replicates --circos-bin-width 10000 --circos-show-chroms I --logging-level INFO --log-file /data/rpbp-predictions-summary.log


Due to the size of the data, we reduce the bin width for the `Circos <http://circos.ca/>`_ plot. We also need to specify which sequences or chromosomes we want to include (by default, only numbered chromosomes and X/x, Y/y are shown). You now have a new directory *analysis* with *profile_construction* and *rpbp_predictions* output.

Launch any of the web applications with

.. code-block:: bash

    docker run -p 8050:8050 --volume `pwd`:/data quay.io/biocontainers/rpbp:3.0.1--py310h30d9df9_0 rpbp-profile-construction-dashboard -c /data/c-elegans-test.yaml --host="0.0.0.0"

.. code-block:: bash

    singularity run --bind `pwd`:/data rpbp.sif rpbp-profile-construction-dashboard -c /data/c-elegans-test.yaml --host="0.0.0.0"

or

.. code-block:: bash

    docker run -p 8050:8050 --volume `pwd`:/data quay.io/biocontainers/rpbp:3.0.1--py310h30d9df9_0 rpbp-predictions-dashboard -c /data/c-elegans-test.yaml --host="0.0.0.0"

.. code-block:: bash

    singularity run --bind `pwd`:/data rpbp.sif rpbp-predictions-dashboard -c /data/c-elegans-test.yaml --host="0.0.0.0"


You then have to open a browser page at the correct address, *e.g.* you see ``Running on http://127.0.0.1:8050``, click on this link, or open a browser page at this address. To navigate the apps is easy, just follow the "hints". Most items are interactive. Press ``CTRL+C`` to quit.
