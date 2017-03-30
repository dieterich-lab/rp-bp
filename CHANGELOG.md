# Change Log
All notable changes to Rp-Bp will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/), 
and this project adheres to [Semantic Versioning](http://semver.org/).

## [in progress]
### Updated
- Predictions report formatting
- Guessing read filtering plot tick marks. See 
  [Issue #60](https://github.com/dieterich-lab/rp-bp/issues/63) for details.

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
