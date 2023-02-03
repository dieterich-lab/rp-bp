Frequently asked questions
==========================


I don't want/I can't install it on my computer. Can I still use **Rp-Bp** without installing the package?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Yes, you can use a Docker or a Singularity container. Simply pull, and you're done! See :ref:`installation_full` for instructions.
Example calls are also given in the user guide and tutorials.

How do I launch the app remotely?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, the application is opened in a browser page on *localhost:8050*. You don't have to actually worry about this. But you can also specify a ``--host`` and a ``--port`` when calling the app, enabling you to launch it from a remote server. In the latter case, however, you have to open a browser page at the correct address. For example, if you use ``--host 123.123.123.123``, then open a page on *http://123.123.123.123:8050/*.

Why is the ORF label not consistent with the host transcript?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In some cases, the ORF label may not be consistent with the host transcript, as reported by the ORF "id" (*transcript_seqname:start-end:strand*). To resolve such seemingly incoherent assignments, compatible transcripts are reported for each ORF in *<genome_name>.orfs-labels.annotated[.orf_note].tab.gz* and shown in the prediction dashboard, see :ref:`apps`.

I have my own alignments, can I use **Rp-Bp**?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The short answer is yes. The pipeline is designed to handle all steps from raw FASTQ files up to the final list of translated Ribo-seq ORFs, but you can start the pipeline from any step. Check the tutorial :ref:`existing-alignment`.
