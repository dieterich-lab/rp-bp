.. _existing-alignment:

How to use existing alignment files
===================================

For this tutorial, we use the *c-elegans-chrI-example* with the *c-elegans.alignments-only.yaml* configuration file.

The **Rp-Bp** pipeline is designed to handle all steps from raw FASTQ files up to the final list of translated Ribo-seq ORFs, but you can start the pipeline from any step.


.. caution::

    Do not use the ``--overwrite`` flag!


For example, the trimming, filtering and alignment steps (performed with Flexbar, Bowtie2, and STAR) could be handled using a different preprocessing strategy and/or different software. In this case, the configuration file must be created as usual, but the ``riboseq_samples`` dictionary only needs to contain the *key* for the samples (the actual path can be left blank or filled with any arbitrary value). This *key* must match the name of any existing files, and these files must be placed at the appropriate location according to the naming convention and nomenclature of **Rp-Bp**.

In our example, the files should be placed (or symlinked) to the following locations:

* Trimmed and quality filtered reads
    *<riboseq_data>/without-adapters/<sample_name>[.note].fastq.gz*

* Reads not aligning to ribosomal sequences
    *<riboseq_data>/without-rrna/<sample_name>[.note].fastq.gz*

* Aligned reads
    *<sample_name>[.note].Aligned.sortedByCoord.out.bam* or
    *<sample_name>[.note].bam* or
    *<sample_name>[.note]-unique.bam* (unless the ``keep_riboseq_multimappers`` configuration option is given)

See the :ref:`rpbp_usage` for more details about required input and output files. Only the "last files" must be in the expected location. For example, if trimming, filtering and alignment has been performed, only the alignment files must be present.


Example
-------

This example shows how to run **Rp-Bp** starting with the alignment files created in the first tutorial :ref:`tutorial_cel`.

.. important::

    You first need to go through the first tutorial, and "run the example dataset" to be able to reproduce this example, as it uses existing alignment files!

Navigate to the same *c-elegans-chrI-example* directory (where you extracted the data).
Create a new directory called *alignments-only*, and symlink existing alignment files

.. code-block:: bash

    mkdir alignments-only && cd alignments-only && mkdir without-rrna-mapping && cd without-rrna-mapping
    # you are now inside */c-elegans-chrI-example/alignments-only/without-rrna-mapping
    # symlink files using absolute or relative path
    # e.g. ln -s <original-example>/without-rrna-mapping/c-elegans-rep-1Aligned.sortedByCoord.out.bam .
    # here we use relative paths
    ln -s ../../without-rrna-mapping/c-elegans-rep-1Aligned.sortedByCoord.out.bam .
    ln -s ../../without-rrna-mapping/c-elegans-rep-2Aligned.sortedByCoord.out.bam .
    # return to c-elegans-chrI-example
    cd ../../


Change the paths in the configuration file, as explained in the first tutorial, or simply modify in place

.. code-block:: bash

    sed -i 's|/path/to/your/c-elegans-example|'`pwd`'|g' c-elegans.alignments-only.yaml

It is now important to also adjust the configuration file by adding the *alignments-only* directory at the end of the path given by the ``riboseq_data`` key in *c-elegans.alignments-only.yaml*

.. code-block:: yaml

    # e.g. if your path is /home/user/data
    riboseq_data: /home/user/data/c-elegans-chrI-example/alignments-only


.. hint::

    Instead of the above commands, you can also create directories using your menu/shortcuts, *etc.*, copy files by hand, and/or use a text editor to modify the configuration file.


Run the pipeline

.. code-block:: bash

    run-all-rpbp-instances c-elegans.alignments-only.yaml --merge-replicates --run-replicates --num-cpus 4 --logging-level INFO --log-file alignments-only.txt

Logging reports a few ``WARNING ... Some input files are missing. Skipping call...``, but the alignment files are eventually found, and the pipeline proceeds from there. You should end up with the same final output as that obtained in the first tutorial.
