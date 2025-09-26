.. _defaults:

Default parameters
==================

The parameters and options decribed below are all optional. All parameters and options have default values that do not normally need to be modified.

.. important::

    **Rp-Bp** parameters are changed via the configuration file, but options for programs such as Flexbar and STAR are handled via command line arguments.
    You do not need to include **Rp-Bp** parameters in the configuration file, unless you wish to change their values.


Flexbar and STAR options
------------------------

Default options are overridden via command line using ``--flexbar-options`` or ``--star-options``. Currently, no options can be passed to Bowtie2.

Flexbar
^^^^^^^

* ``max-uncalled`` Default: 1.
* ``pre-trim-left`` Default: 0.
* ``qtrim-format`` Default: sanger.
* ``qtrim`` Default: TAIL.
* ``qtrim-threshold`` Default: 10.
* ``zip-output`` Default: GZ.

STAR
^^^^

* ``readFilesCommand`` Default: zcat (gunzip -c for macOS).
* ``limitBAMsortRAM`` Default: 0 (set to ``--mem`` at run-time).
* ``alignIntronMin`` Default: 20.
* ``alignIntronMax`` Default: 100000.
* ``outFilterMismatchNmax`` Default: 1.
* ``outFilterMismatchNoverLmax`` Default: 0.04.
* ``outFilterType`` Default: BySJout.
* ``outFilterIntronMotifs`` Default: RemoveNoncanonicalUnannotated.
* ``outSAMattributes`` Default: AS NH HI nM MD.
* ``outSAMtype`` Default: BAM SortedByCoordinate.
* ``sjdbOverhang`` Default: 33.
* ``seedSearchStartLmaxOverLread`` Default: 0.5.
* ``winAnchorMultimapNmax`` Default: 100.


Rp-Bp parameters
----------------

* ``keep_riboseq_multimappers`` If ``True`` in the configuration file, then multimapping riboseq reads *will not* be removed. They will be treated as "normal" reads in every place they map, *i.e.* the weight of the read will not be distributed fractionally.
* ``models_base`` The path to the compiled models, if installed in a different location. The models are included with the source distribution and compiled as part of the installation. *Do not change this, unless you know what you are doing!*


Shared MCMC parameters
^^^^^^^^^^^^^^^^^^^^^^

* ``seed`` The random seed for the MCMC sampling, used for periodicity estimation and translation prediction. Default: 8675309.
* ``chains`` The number of chains to use in the MCMC sampling, used for periodicity estimation and translation prediction. Default: 2


Metagene and periodicity estimation parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*  ``metagene_start_upstream`` The number of bases upstream of the translation initiation site to begin constructing the metagene profile. Default: 300.
*  ``metagene_start_downstream`` The number of bases downstream of the translation initiation site to end the metagene profile. Default: 300.
*  ``metagene_end_upstream`` The number of bases upstream of the translation termination site to begin constructing the metagene profile. Default: 300.
*  ``metagene_end_downstream`` The number of bases downstream of the translation termination site to end the metagene profile. Default: 300.
*  ``periodic_offset_start`` The position, relative to the translation initiation site, to begin calculating periodicity Bayes factors. Default: -20 (inclusive).
*  ``periodic_offset_end`` The position, relative to the translation initiation site, to stop calculating periodicity Bayes factors. Default: 0 (inclusive).
*  ``metagene_profile_length`` The length of the profile to use in the models. ``metagene_profile_length`` + ``periodic_offset_end`` must be consistent with the length of the extracted metagene profile. Default: 21.
*  ``metagene_iterations`` The number of iterations to use for each chain in the MCMC sampling. Default: 500 (includes warmup).
*  ``min_metagene_profile_count`` Read lengths with fewer than this number of reads will not be used. Default: 1000.
*  ``min_metagene_bf_mean`` If ``max_metagene_bf_var`` and ``min_metagene_bf_likelihood`` are None (null in YAML), this is taken as a hard threshold on the estimated Bayes factor mean. Default: 5.
*  ``max_metagene_bf_var`` A hard threshold on the estimated Bayes factor variance. Default: None.
*  ``min_metagene_bf_likelihood`` A threshold on the likelihood of periodicity. Default: 0.5.


.. note::

    A profile is periodic if [P(bf > ``min_metagene_bf_mean``)] > ``min_metagene_bf_likelihood``. By default, we do not filter on the variance. If given, then both filters are applied and the result is the intersection.


Fixed lengths and offsets
^^^^^^^^^^^^^^^^^^^^^^^^^

* ``use_fixed_lengths`` If ``True`` in the configuration file, fixed values given by ``lengths`` and ``offsets`` are used (no periodicity estimation).
* ``lengths`` A list of read lengths to use for creating the profiles if the ``use_fixed_lengths`` option is ``True``. Presumably, these are lengths that have periodic metagene profiles.
* ``offsets``  The P-site offset to use for each read length specifed by ``lengths`` if the ``use_fixed_lengths`` option is ``True``. The number of offsets must match the number of lengths, and they are assumed to match. For example ``lengths``:  [26, 29] with ``offsets``: [9, 12] means only reads of lengths 26 bp and 29 bp are used to create the profiles. The 26 bp reads will be shifted by 9 bp in the 5' direction, while reads of length 29 bp will be shifted by 12 bp.


Smoothing parameters
^^^^^^^^^^^^^^^^^^^^

* ``smoothing_fraction`` The fraction of the data used when estimating each y-value for LOWESS. Default: 0.2.
* ``smoothing_reweighting_iterations`` The number of residual-based reweightings to perform for LOWESS. See the `statsmodels documentation <https://www.statsmodels.org>`_. Default: 0.


Translation prediction parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``orf_min_length_pre`` ORFs with length < ``orf_min_length_pre`` (nucleotides) are not processed. Default: 0 (ignore option).
* ``orf_max_length_pre`` ORFs with length > ``orf_max_length_pre`` (nucleotides) are not processed. Default: 0 (ignore option).
* ``orf_min_length`` Only ORFs with length > ``orf_min_length`` (nucleotides) are kept in the final set. Default: 8.
* ``orf_min_profile_count_pre`` ORF with profile sum < ``orf_min_profile_count_pre`` are not processed. Default: 5.
* ``orf_min_profile_count`` Only ORFs with profile sum > ``orf_min_profile_count`` are kept in the final set. Default: None.
* ``translation_iterations`` The number of iterations to use for each chain in the MCMC sampling. Default: 500 (includes warmup).
* ``min_bf_mean`` If ``max_bf_var`` and ``min_bf_likelihood`` are None (null in YAML), this is taken as a hard threshold on the estimated Bayes factor mean. Default: 5.
* ``max_bf_var`` A hard threshold on the estimated Bayes factor variance. Default: None.
* ``min_bf_likelihood`` A threshold on the likelihood to select an ORF as translated. Default: 0.5.
* ``chisq_alpha`` For the chi-square test, this value is first Bonferroni corrected based on the number of ORFs which pass the smoothing filters. It is then used as the significance threshold to select translated ORFs. Default: 0.01.


.. note::

    A Ribo-seq ORF is translated if [P(bf > ``min_bf_mean``)] > ``min_bf_likelihood``. By default, we do not filter on the variance. If given, then both filters are applied and the result is the intersection.


.. attention::

    Chi-square values are reported, but they are not used for prediction, unless the ``chi_square_only`` flag is present in the configuration file, in which case the translation models are not fit to the data, and the posterior distributions are not estimated. This is mostly kept for historical reasons, and may eventually be removed.
