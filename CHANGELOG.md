# Change Log
All notable changes to Rp-Bp will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/), 
and this project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased] - started 2019-10-16

## [2.0.0] 2019-05-24

This is a major version upgrade due to changes in API and package dependencies. 

### Changed
- Configure setup() using setup.cfg, add `pbio`, update install via requirements.

### Added
- Config file for RTD
- Module with all defaults parameters, handling options for external dependencies
- Add option to force the pickling of the Stan models during installation 
    [ISSUE #88](https://github.com/dieterich-lab/rp-bp/issues/88).

### Fixed
- Build RTD documentation (add `autodoc_mock_imports`, empty `docs/requirements.txt`)
- Missing parameters in function calls, hard coded defaults
- [YAMLLoadWarning](https://msg.pyyaml.org/load)
- Fall back to multiprocessing in `estimate_orf_bayes_factors`, and 
    `extract_orf_profiles`
- Stan assignment operator in models
- ORF labels in `prepare_rpbp_genome`, `label_orfs` and `extract_orf_coordinates`, see
    [ISSUE #96](https://github.com/dieterich-lab/rp-bp/issues/96).
    
### Removed
- Transcript exons file (unused)
 
    
## [1.1.12] - 2018-03-14
### Changed
- Example dataset updated and moved to new address. Updated documentation.
- `predict_translated_orfs` and `select_final_prediction_set` modified such 
    that Rp-chi pipeline is used *only* if `chi_square_only` is in config file,
    otherwise Rp-Bp is used.
- Handling of options passed to flexbar via arguments to `run_all_rpbp_instances`,
    `run_rpbp_pipeline`, `create_orf_profiles` and `create_base_genome_profile`. 
- Handling of additional options passed to star via arguments, mainly 
    to `create_base_genome_profile`. 

### Removed
- Removed config option `seqids_to_keep` for compatibility with `extract_metagene_profiles`
    that does not use this option anymore. This has also been removed from the example. 
    Updated documentation.
- Removed fields `delta_l` and `delta_h` from output in `estimate_orf_bayes_factors`.
- Removed deprecated function calls from 
    `analysis/profile_construction/get_all_read_filtering_counts`,
    `prepare_rpbp_genome`, `create_base_genome_profile`

### Added
- Utility for supressing pystan (or other compiled function) output. For more
    details see [Issue #10](https://github.com/dieterich-lab/rp-bp/issues/10).
- Added *limited support* for GFF3 specifications, in particular this requires removing the
    STOP codons from the CDSs and passing the right option to `STAR`. 
    `Prepare_rpbp_genome` and `create_base_genome_profile` modified.
- If using de_novo, all annotations are concatenated automatically in preparation of
    running the main pipeline.

### Fixed
- Minor changes to plotting options, typos and/or redundant features in 
    `rpbp/analysis/profile_construction/create_rpbp_preprocessing_report` and 
    `rpbp/analysis/profile_construction/visualize_metagene_profile_bayes_factor`. 
    See [ISSUE #87](https://github.com/dieterich-lab/rp-bp/issues/87). In progress.
- Added default *sjdbOverhang* option to star during mapping.

## [1.1.11] - 2017-12-08
### Removed
- Removed deprecated function calls from `analysis/rpbp_predictions/add_mygene_info_to_orfs`, 
    `analysis/find_differential_micropeptides`, `analysis/proteomics/get_orf_peptide_matches` 
    and `analysis/proteomics/filter_nonunique_peptide_matches`

### Fixed
- In `create_rpbp_preprocessing_report.create_figures` call to 
    `filenames.get_riboseq_profiles` modified to fix file name reference issue.
- In `create_rpbp_preprocessing_report` some figures were referenced before they were created,
    as a *temporary fix*, the function `create_figures` is now called earlier.
- Minor typos in some of the analysis scripts.

### Added
- Flag added to differentiate between *exons* file (list of exons for any given transcript),
    and *orfs* file (technically, a list of exons per orf). 
    See [Issue #59](https://github.com/dieterich-lab/rp-bp/issues/59) for more details.

### Updated
- Dependencies to new repo locations.

## [1.1.10] - 2017-10-27
### Added
- Reference to helper script to install prerequisite programs. See 
  [Issue #75](https://github.com/dieterich-lab/rp-bp/issues/75) for more
  details.

### Updated
- Dependencies to new `pymisc-utils` and `pybio-utils` repos. See
  [Issue #74](https://github.com/dieterich-lab/rp-bp/issues/74) for more
  details.

- Version specification of prereqs

### Removed
- `cluster-subcodon-counts`. See [Issue #78](https://github.com/dieterich-lab/rp-bp/issues/78)
  for more details.

## [1.1.9] - 2017-06-15
### Added
- Documentation for using custom alignment files. Please see
  [Issue #73](https://github.com/dieterich-lab/rp-bp/issues/73) for more
  details.

### Fixed
- Documentation links to old ipython notebooks. Please see
  [Issue #71](https://github.com/dieterich-lab/rp-bp/issues/71) for more
  details.

## [1.1.8] - 2017-05-22
### Fixed
- Incorrect read filtering report. See [Issue #69](https://github.com/dieterich-lab/rp-bp/issues/69)
  for more details.

- Broken download link for example.

## [1.1.7] - 2017-04-03
### Fixed
- Incorrect alignment of ORF data frame and profile matrix. See
  [Issue #54](https://github.com/dieterich-lab/rp-bp/issues/54) for more details.
- ORF coordinates for start codon on exon boundary. See
  [Issue #64](https://github.com/dieterich-lab/rp-bp/issues/64) after the
  reopening for more details.


## [1.1.6] - 2017-03-31
### Updated
- Predictions report formatting
- Guessing read filtering plot tick marks. See 
  [Issue #60](https://github.com/dieterich-lab/rp-bp/issues/60) for details.

### Fixed
- ORF coordinates when start and stop codons fall on an exon boundary. See
  [Issue #64](https://github.com/dieterich-lab/rp-bp/issues/64) for more details.

## [1.1.5] - 2017-03-29
### Added
- Documentation for the analysis scripts
- Documentation for creating reference indices

### Fixed
- Checks that annotation files exist before beginning index creation. See 
  [Issue #62](https://github.com/dieterich-lab/rp-bp/issues/62)
- Remove unnecessary `--tmp` argument for pipeline call to
  `predict-translated-orfs`. See [Issue #63](https://github.com/dieterich-lab/rp-bp/issues/63)
  for more details.

## [1.1.4] - 2017-02-22
### Removed
- Unnecessary [analysis] package. Depending on how the package was installed,
  this sometimes resulted in the scripts not being found. See [Issue #56]
  (https://github.com/dieterich-lab/rp-bp/issues/56) for details.

## [1.1.3] - 2017-02-21
### Fixed
- Handling of Stan model paths which include spaces. See [Issue #36](https://github.com/dieterich-lab/rp-bp/issues/36) 
  for more details. This bug was reintroduced with the commit:
  "UPD resolve conflicts from dev for version 1.1"

## [1.1.2] - 2017-02-15
### Fixed
- Backward-compatiblity syntax in `estimate-orf-bayes-factors`. 
  See [Issue #51](https://github.com/dieterich-lab/rp-bp/issues/51#issuecomment-280024994)
  for more details.

### Added
- `__version__` and `__version_info__` for `rpbp` package

## [1.1.1] - 2017-02-14
### Updated
- Documentation to use pandoc for converting markdown to html

### Fixed
- Broken download link to example files

## [1.1.0] - 2017-02-14

This is a rather significant update to the codebase and includes many small 
changes. Primarily, though, it addresses many installation and system issues.

In particular, thanks to Tonu Margus ([@tmargus](https://github.com/tmargus))
for working through many of these issues.

### Added
- Installation instructions for anaconda

### Fixed
- Various installation issues related to calling ``pip3`` from ``setup.py``.
  Installation of private repositories is now handled via a ``requirements.txt``
  file.

### Removed
- All uses of (py)bedtools

## [1.0.1] - 2017-02-01
### Fixed
- Handling of Stan model paths which include spaces. In particular, this was a
  problem on OSX. See [Issue #36](https://github.com/dieterich-lab/rp-bp/issues/36)
  for more details.


## [1.0.0] - 2016-09-05
The initial version, which implements everything in the paper. It was only used
on debian-based systems.
