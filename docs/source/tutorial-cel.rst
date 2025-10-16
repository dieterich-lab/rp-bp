.. _tutorial_cel:

Run the ORF discovery pipeline
==============================

For this tutorial, we use the *c-elegans-chrI-example* with the *c-elegans-test.yaml* configuration file.

How to run the pipeline
-----------------------

To run the pipeline, change the paths in the configuration file to point to the location where you extracted the data, *e.g.* */home/user/data/c-elegans-chrI-example*. You can do this using a text editor, or simply by modifying the file in place

.. code-block:: bash

    sed -i 's|/path/to/your/c-elegans-example|'`pwd`'|g' c-elegans-test.yaml


.. important::

    Default parameters were modified for the example and included in the configuration file. If you use this configuration file as a general template for your data, do not forget to remove everything below the line "REMOVE BELOW THIS LINE IF YOU USE THIS CONFIGURATION FILE AS TEMPLATE FOR YOUR DATA".


You can now create the genome indices and annotations. For this small example, it is important to downscale the STAR ``--genomeSAindexNbases`` parameter as follows


.. code-block:: bash

    prepare-rpbp-genome c-elegans-test.yaml --star-options "--genomeSAindexNbases 10" --num-cpus 4 --logging-level INFO --log-file rpbp-genome.log

The file *rpbp-genome.log* contains logging output for the reference preprocessing. You now have a new directory called *WBcel235.79.chrI* with genome indices and annotations.

Finally, run the ORF discovery pipeline


.. code-block:: bash

    run-all-rpbp-instances c-elegans-test.yaml --merge-replicates --run-replicates --keep-intermediate-files --num-cpus 4 --logging-level INFO --log-file rpbp-pipeline.log

The file *rpbp-pipeline.log* contains logging output for the different processing steps. You now have four new directories (*with-*, *without-*) including output from Flexbar, Bowtie2, and STAR, and directories with **Rp-Bp** output: *metagene-profiles*, *orf-profiles*, and *orf-predictions*. The *orf-predictions* include the output for each sample *c-elegans-rep-1* and *c-elegans-rep-2* as well as for the merged replicates *c-elegans-test*.

How to summarize the results and launch the apps
------------------------------------------------

Prepare the summary output for the *profile construction dashboard*

.. code-block:: bash

    summarize-rpbp-profile-construction c-elegans-test.yaml --num-cpus 4 --logging-level INFO --log-file rpbp-profile-summary.log


and for the *predictions dashboard*

.. code-block:: bash

    summarize-rpbp-predictions c-elegans-test.yaml --no-replicates --circos-bin-width 10000 --circos-show-chroms I --logging-level INFO --log-file rpbp-predictions-summary.log

Due to the size of the data, we reduce the bin width for the `Circos <http://circos.ca/>`_ plot. We also need to specify which sequences or chromosomes we want to include (by default, only numbered chromosomes and X/x, Y/y are shown). You now have a new directory *analysis* with *profile_construction* and *rpbp_predictions* output.

Launch any of the web applications with

.. code-block:: bash

    rpbp-profile-construction-dashboard -c c-elegans-test.yaml

or

.. code-block:: bash

    rpbp-predictions-dashboard -c c-elegans-test.yaml


To navigate the apps is easy, just follow the "hints". Most items are interactive. Press ``CTRL+C`` to quit. Check :ref:`apps` for more details.

.. attention::

    For the apps only, the configuration file is passed using a (required) named argument ``-c/--config CONFIG``.
