.. _running_rpbp:

How to run the pipeline
=======================

See `ORF profile construction`_ and `Translation prediction`_ for a short description of required input and expected output. **Rp-Bp** output files are written in the `BED <https://www.ensembl.org/info/website/upload/bed.html>`_, `FASTA <https://en.wikipedia.org/wiki/FASTA_format>`_, `sparse matrix market (MTX) <http://math.nist.gov/MatrixMarket/formats.html>`_, or CSV format. Output from Flexbar, Bowtie2, and STAR are written in FASTQ or `BAM <https://samtools.github.io/hts-specs/>`_ formats.

.. attention::

    All Ribo-seq samples (including biological replicates) in the configuration file must be from the same organism and use the same ``genome_base_path``, ``star_index``, ``ribosomal_index``, *etc.* Samples from different organisms or using different annotations must be "split" into different configuration files, and run separately.

.. _rpbp_usage:

General usage
-------------

.. code-block:: bash

    # Only create the ORF profiles (estimate periodicity) for QC.
    run-all-rpbp-instances --profiles-only [options] config

    # Run the ORF discovery pipeline for all samples in the configuration file (only samples, i.e. do not merge the replicates).
    run-all-rpbp-instances [options] config

    # Run the ORF discovery pipeline for all samples in the configuration file, merge the replicates, and make predictions for merged replicates.
    run-all-rpbp-instances --merge-replicates --run-replicates [options] config

For all options, consult the API for :ref:`api_prepare`. See also :ref:`howto_config`.

.. tip::

    To perform read filtering quality control (QC), use the ``-k/--keep-intermediate-files`` option, see :ref:`apps`.

----

.. _biological_replicates:

More about biological replicates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To merge biological replicates, use the ``--merge-replicates`` flag. This tells **Rp-Bp** to use the combined ORF profiles to estimate the Bayes factors. If the ``--merge-replicates`` flag is given, predictions will not be made for the individual samples, unless the ``--run-replicates`` flag is also given, in which case predictions are made for both the merged replicates as well as the individual samples. This is how you would prepare a configuration file for four samples of two different conditions:

.. code-block:: yaml

    # example of 4 samples: 2 controls and 2 conditions
    # <sample_name> below is replaced by ctrl1, ctrl2, cond1, and cond2
    # <condition_name> below is replaced by ctrl and cond

    riboseq_samples:
     ctrl1: /path/to/sample1.fastq.gz
     ctrl2: /path/to/sample2.fastq.gz
     cond1: /path/to/sample3.fastq.gz
     cond2: /path/to/sample4.fastq.gz

    riboseq_biological_replicates:
     ctrl:
      - ctrl1
      - ctrl2
     cond:
      - cond1
      - cond2

    # fancy names to use for downstream analyses
    riboseq_sample_name_map:
     ctrl1: Ctrl-1
     ctrl2: Ctrl-2
     cond1: Cond-1
     cond2: Cond-2

    riboseq_condition_name_map:
     ctrl: Ctrl
     cond: Cond


.. _run_profile:

ORF profile construction
------------------------

To run the periodicity estimation only, pass the ``--profiles-only`` option.

.. note::

    This part of the pipeline uses Flexbar, Bowtie2, and STAR to process and align Ribo-seq reads, however you can estimate periodicity (and predict translation events) using your own existing alignment files in BAM format, see :ref:`existing_alignment`.


Required input
^^^^^^^^^^^^^^

All the input files are given in the configuration file.

Expected output
^^^^^^^^^^^^^^^

The base path for the following files is: *<riboseq_data>/without-adapters*

* *<sample_name>[.note].fastq.gz*. Output from Flexbar (adapters and low-quality reads removed).

The base path for the following files is: *<riboseq_data>/with-rrna*

* *<sample_name>[.note].fastq.gz*. Output from Bowtie2 (reads aligning to the ribosomal index, kept for quality control, but unused).

The base path for the following files is: *<riboseq_data>/without-rrna*

* *<sample_name>[.note].fastq.gz*. Output from Bowtie2 (after *in-silico* rRNA removal, used for genome alignment).

The base path for the following files is: *<riboseq_data>/without-rrna-mapping*

* *<sample_name>[.note].bam*. A sorted BAM file with genome alignments (the *Aligned.sortedByCoord.out.bam* STAR output).
* *<sample_name>[.note]-unique.bam*. A sorted BAM file with unique alignments (multimapping reads removed).


.. note::

    If ``keep_riboseq_multimappers`` is ``True`` in the configuration file, then there will be no *-unique* files. In general, we do not recommend to keep multimappers.


The base path for the following files is: *<riboseq_data>/metagene-profiles*

* *<sample_name>[.note][-unique].metagene-profile.csv.gz*. A CSV file with the metagene profiles constructed from aligned reads (given by the "position" or offset and "count" columns) for all read lengths ("length" column). It include profiles for the annotated translation initiation sites and translation termination sites ("type" column).
* *<sample_name>[.note][-unique].metagene-periodicity-bayes-factors.csv.gz*. A CSV file with the model outputs and Bayes factor estimates for all P-site offsets and read lengths.
* *<sample_name>[.note][-unique].periodic-offsets.csv.gz*. A CSV file with the best P-site offset for each read length. All read lengths are included, even if the estimates do not meet the prediction criteria (filtering occurs on the fly).

The base path for the following files is: *<riboseq_data>/orf-profiles*

* *<sample_name>[.note][-unique].length-<lengths>.offset-<offsets>.profiles.mtx.gz*. A MTX file with the profiles for all putative Ribo-seq ORFs ("orf_num", "orf_position", *i.e.* position within the ORF, and "read_count"). The matrix market format uses base-1 indexing!
* *<condition_name>[.note][-unique].profiles.mtx.gz*. Same as above for condition, if using ``--merge-replicates``.


.. _run_predictions:

Translation prediction
----------------------

Without the ``--profiles-only`` option, the pipeline will predict which ORFs show evidence of translation, using only the periodic footprint lengths. The ``--merge-replicates`` options is used to predict translation events in merged profiles.

.. tip::

    If you first created profiles and estimated periodicity using the ``--profiles-only`` option, you can decide to continue with the translation prediction step at a later stage. You only have to ``run-all-rpbp-instances <config> [--merge-replicates] [--run-replicates]`` using the same configuration file. Steps for which output files already exists will be skipped, unless the ``--overwrite`` option is set.


Required input
^^^^^^^^^^^^^^

In addition to input files given in the configuration file, metagene and ORF profile output files are required (see above). If the pipeline is run sequentially, you do not normally have to worry about the intermediate output.


Expected output
^^^^^^^^^^^^^^^

The base path for the following files is: *<riboseq_data>/orf-predictions*

* *<sample_name>[.note][-unique].length-<lengths>.offset-<offsets>.bayes-factors.bed.gz*. A BED12+ file with model outputs for all Ribo-seq ORFs. Additional columns include the ORF number, ORF length, model outputs, Bayes factor mean and variance, and P-site coverage across 3 frames.
* *<sample_name>[.note][-unique].length-<lengths>.offset-<offsets>[.filtered].predicted-orfs.bed.gz*. Same format as above, with the predicted translation events. **This file contains the translated Ribo-seq ORFs**.
* *<sample_name>[.note][-unique].length-<lengths>.offset-<offsets>[.filtered].predicted-orfs.dna.fa*. A FASTA file with the predicted translation events. The FASTA header matches the "id" column in the corresponding BED file. **This file contains the DNA sequence for each translated Ribo-seq ORF**.
* *<sample_name>[.note][-unique].length-<lengths>.offset-<offsets>[.filtered].predicted-orfs.protein.fa*. A FASTA file with the predicted translation events. The FASTA header matches the "id" column in the corresponding BED file. **This file contains the protein sequence for each translated Ribo-seq ORF**.

* *<condition_name>[.note][-unique].bayes-factors.bed.gz*. Same as above for condition, if using ``--merge-replicates``.
* *<condition_name>[.note][-unique][.filtered].predicted-orfs.bed.gz*. Same as above for condition, if using ``--merge-replicates``.
* *<condition_name>[.note][-unique][.filtered].predicted-orfs.dna.fa*. Same as above for condition, if using ``--merge-replicates``.
* *<condition_name>[.note][-unique][.filtered].predicted-orfs.protein.fa*. Same as above for condition, if using ``--merge-replicates``.

.. attention::

    Translation events are predicted using Bayesian model selection. Our model does not distinguishes between overlapping ORFs. To select the best overlapping ORF among a group of overlapping ORFs, we first select the longest ORF, then the highest Bayes factor. This is referred to as the *filtered* predictions.

    In previous versions, both *filtered* and *unfiltered* (including all overlapping ORFs) predictions were written to file. In general, we recommend to use *filtered* predictions. Unless the ``--write-unfiltered`` option is used, **Rp-Bp** now only outputs the *filtered* predictions. If using ``--write-unfiltered``, *unfiltered* predictions are also written to file, without the *[.filtered]* flag. Hence to avoid confusion with older results, the *filtered* predictions have kept the *[.filtered]* flag.


.. note::

    If *smoothing parameters* (see :ref:`defaults`) are given in the configuration file, the following string *.frac-<smoothing_fraction>.rw-<smoothing_reweighting_iterations>* is also added to the file names. Default values (unless they are explicitly given in the configuration file) are not written.
