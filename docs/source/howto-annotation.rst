.. _howto_annotation:

How to prepare annotations
==========================

Run this once for any given reference genome, as long as the package version and its dependencies remain the same.

.. _genome_usage:

General usage
-------------

.. code-block:: bash

    prepare-rpbp-genome [options] config

For all options, consult the API for :ref:`api_prepare`. See also :ref:`howto_config`.


Required input
--------------

Reference annotations
^^^^^^^^^^^^^^^^^^^^^

A reference GTF annotation (with a gtf extension) for your genome, with *transcript*, *exon*, and *CDS* features (*start_codon* and *stop_codon* features are not required). The attribute field must include *transcript_id* and *gene_id*, and optionally *transcript_biotype*, *gene_name*, and *gene_biotype*. The annotation must match the version of the reference genome sequence.

.. caution::

    The GTF/GFF2 format does not include the stop codon in the terminal exon, *i.e.* the stop codon is not included in the CDS.


Reference genome sequence
^^^^^^^^^^^^^^^^^^^^^^^^^

The input FASTA file for your genome including top level sequences, but excluding haplotypes, *e.g.* the *primary assembly* file from `Ensembl <https://www.ensembl.org/info/data/ftp/index.html>`_. The identifiers must match those in the GTF annotation file.


Ribosomal sequence
^^^^^^^^^^^^^^^^^^

A separate FASTA file with ribosomal DNA (rDNA) sequences which are generally not included in the genome assembly. This is used for *in silico* rRNA removal. This file can also include other sequences to filter out, *.e.g* pre-rRNA, tRNA. We typically include the following

* The large and small ribosomal subunit sequences, *e.g.* from NCBI.
* The genomic tRNA sequences *e.g.* from `GtRNAdb <http://gtrnadb.ucsc.edu>`_.
* Mt_rRNA, Mt_tRNA and rRNA genes from BioMart. Select options for the "Gene type" filter. For "Attributes", select "Sequences", and then specifically "Exon sequences". Additionally, including the "Gene type" in the header can be helpful for quality control.


Expected output
---------------

Output files are written in the `BED <https://www.ensembl.org/info/website/upload/bed.html>`_,  `FASTA <https://en.wikipedia.org/wiki/FASTA_format>`_, or TAB-delimited formats.


* *<ribosomal_index>*. The `Bowtie 2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_ index.

* *<star_index>*. The `STAR <https://github.com/alexdobin/STAR>`_ index.


The base path for the following file is: *<genome_base_path>*


* *<genome_name>.annotated.bed.gz*. Annotated transcripts in BED12+ format.


The base path for the following files is: *<genome_base_path>/transcript-index*


* *<genome_name>.transcripts.annotated.fa*. Annotated transcript sequences in FASTA format.

* *<genome_name>.orfs-genomic[.orf_note].bed.gz*. Putative Ribo-seq ORFs in BED12+ format. See :ref:`rpbp_genome` for details.

* *<genome_name>.orfs-exons[.orf_note].bed.gz*. Putative Ribo-seq ORF exons in BED6+ format. The extra columns are *exon_index*, giving the order of the exon in the transcript, and *transcript_start*, giving the start position of that index in transcript coordinates.

* *<genome_name>.orfs-labels[.orf_note].tab.gz*. Putative Ribo-seq ORF categories and compatible transcripts in TAB-delimited format. See :ref:`rpbp_genome` for details.


.. note::

    If a ``de_novo_gtf`` file is provided, the last three output files are split into *annotated* and *de-novo*. The files used by the pipeline, as described above, are the "concatenation" of the respective *annotated* and *de-novo* files. In addition, *de-novo* transcript BED and FASTA files are created. A new GTF file, created by concatenating ``gtf`` and ``de_novo_gtf``, is written to ``genome_base_path``. This file may contain repeated features. For mapping with STAR, this is generally not a problem, but it can be for abundance estimation, *cf.* `Ribotools documentation <https://ribotools.readthedocs.io/en/latest/howto-run.html#general-usage>`_.
