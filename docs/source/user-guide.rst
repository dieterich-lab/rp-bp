.. _user_guide:

User guide
==========

You want to "quality control" your Ribo-seq samples for downstrean analyses, or run the Ribo-seq ORF discovery pipeline? First, you need to prepare genome indices and annotations for your organism. This has to be done once for any given reference genome and annotation. Consult :ref:`howto_annotation`.

The pipeline itself consists of two "modules": the *ORF profile construction*, where periodic read lengths and ribosome P-site offsets are inferred from the data; and the *translation prediction*, where translation events are predicted. Consult :ref:`running_rpbp`.

.. toctree::
   :maxdepth: 1

   howto-config
   howto-annotation
   rpbp-genome
   howto-run
   howto-qc
   defaults
