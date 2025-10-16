.. _howto_config:

How to add a config file
========================

Create a single YAML configuration file with the following keys:


* ``gtf`` *(required, input)* Path to the reference annotations (format GTF2).
* ``fasta`` *(required, input)* Path to the reference genome sequence (format FASTA).
* ``ribosomal_fasta`` *(required, input)* Path to the ribosomal sequence (format FASTA).

* ``de_novo_gtf`` *(optional, input)* An additional GTF containing annotations constructed from a *de novo* assembly. See :ref:`denovo`.
* ``start_codons`` *(optional, input)* A list of strings to use as start codons. Default: ``ATG``.
* ``stop_codons`` *(optional, input)* A list of strings to use as stop codons. Default: ``TAA``, ``TGA``, ``TAG``.

* ``genome_name`` *(required, output)* A descriptive name to use for the created files.
* ``genome_base_path`` *(required, output)* Output path (directory) for the transcript fasta and ORFs.
* ``ribosomal_index`` *(required, output)* Output path (directory/filename) for the `Bowtie 2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_ index.
* ``star_index`` *(required, output)* Output path (directory) for the `STAR <https://github.com/alexdobin/STAR>`_ index.
* ``orf_note`` *(optional, output)* An additional description used in the filenames. It should not contain spaces or special characters.


* ``riboseq_samples`` *(required, input)* A dictionary *key: value*, where *key* is used to construct filenames, and *value* is the full path to the FASTQ.gz file for a given sample. The *key* should not contain spaces or special characters.

* ``riboseq_biological_replicates`` *(optional, input)* A dictionary *key: value*, where *key* is a condition, and *value* contains all samples which are replicates of the condition. Items of the *value* list must match the ``riboseq_samples`` *key*.
* ``adapter_file`` *(optional, input)* Path to adapter sequences (FASTA file) to be removed by `Flexbar <https://github.com/seqan/flexbar>`_.
* ``adapter_sequence`` *(optional, input)* A single adapter sequence to be removed. If both ``adapter_file`` and ``adapter_sequence`` are given, the former has precedence.

* ``riboseq_data`` *(required, output)* The base output location for all created files.

* ``riboseq_sample_name_map`` *(optional, output)* A dictionary *key: value*, where *key* is the same as ``riboseq_samples`` *key*, and *value* is a fancy name for *key* to use in downstream analyses.
* ``riboseq_condition_name_map`` *(optional, output)* A dictionary *key: value*, where *key* is the same as ``riboseq_biological_replicates`` *key*, and *value* is a fancy name for *key* to use in downstream analyses.
* ``project_name`` *(optional, output)* An additional description used in the filenames for downstream analyses. It should not contain spaces or special characters. See :ref:`apps`.
* ``note`` *(optional, output)* An additional description used in the filenames. It should not contain spaces or special characters.


To download an example configuration file, check the test Ribo-seq dataset included with the :ref:`all_tutorials`. To change the default parameters, see :ref:`defaults`. See also :ref:`biological_replicates` for an example of how to use biological replicates.
