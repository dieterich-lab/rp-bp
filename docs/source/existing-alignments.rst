.. _existing-alignment:

How to use existing alignment files
===================================

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
    *<sample_name>[.note].bam*
    *<sample_name>[.note]-unique.bam* (unless the ``keep_riboseq_multimappers`` configuration option is given)
    
 
See the :ref:`rpbp_usage` for more details about required input and output files. Only the last files must be in the expected location. For example, if trimming, filtering and alignment has been performed, only the alignment files must be present. The pipeline will issue warning messages that some files are missing, but will continue once it finds the alignment files.


Example
-------

This example shows how to run **Rp-Bp** starting with the alignment files for the `Tutorial <tutorial.html>`_, using the *c-elegans.alignments-only.yaml*.

.. important::

    You first need to go through the tutorial, and "run the example dataset" to be able to reproduce this example, as it uses existing alignment files!
    You also need to update the paths in the *c-elegans.alignments-only.yaml* configuration file!
    

First create the necessary directory structure under *<riboseq_data>*, and symlink existing files 

.. code-block:: bash
    
    # we are inside <riboseq_data> given by c-elegans.alignments-only.yaml
    mkdir without-rrna-mapping && cd without-rrna-mapping
    
    # <original-example> is the base path <riboseq_data> from c-elegans-test.yaml, i.e.
    # the output from the tutorial
    ln -s <original-example>/without-rrna-mapping/c-elegans-rep-1.bam c-elegans-rep-1.bam
    ln -s <original-example>/without-rrna-mapping/c-elegans-rep-1.bam.bai c-elegans-rep-1.bam.bai
    ln -s <original-example>/without-rrna-mapping/c-elegans-rep-2.bam c-elegans-rep-2.bam
    ln -s <original-example>/without-rrna-mapping/c-elegans-rep-2.bam.bai c-elegans-rep-2.bam.bai


Run the pipeline

.. code-block:: bash

    run-all-rpbp-instances c-elegans.alignments-only.yaml --merge-replicates --run-replicates --num-cpus 4 --logging-level INFO --log-file alignments-only.txt
    

Logging reports various missing files, *e.g.*

.. code-block:: bash

    WARNING  misc.shell_utils 2017-06-15 19:02:35,224 : Some input files ['/this/does/not/matter.rep-1.fastq.gz'] are missing. Skipping call: flexbar --qtrim-format sanger --max-uncalled 1  -n 2 --zip-output GZ -r /this/does/not/matter.rep-1.fastq.gz -t /home/bmalone/python-projects/rp-bp/c-elegans.alignments-only/without-adapters/c-elegans-rep-1 --pre-trim-left 0


but the alignment files are eventually found, and the pipeline proceeds from there. 

