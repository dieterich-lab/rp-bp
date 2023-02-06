# Ribosome profiling with Bayesian predictions (Rp-Bp)

Ribosome profiling (Ribo-seq) is an RNA-sequencing-based readout of RNA translation. Isolation and deep-sequencing of ribosome-protected RNA fragments (ribosome footprints) provides a genome-wide snapshot of the translatome at sub-codon resolution. **Rp-Bp** is an unsupervised Bayesian approach to predict translated open reading frames (ORFs) from ribosome profiles. **Rp-Bp** can be used for ORF discovery, or simply to estimate periodicity in a set of Ribo-seq samples. When used for ORF discovery, **Rp-Bp** automatically classifies ORFs into different biotypes or categories, relative to their host transcript.

**Rp-Bp** comes with two _interactive dashboards_ or _web applications_, one for read and periodicity quality control, the other to facilitate Ribo-seq ORFs discovery.

![rpbp](docs/source/_static/logo-rpbp-dark.png)

---

## Documentation

Consult the [user guide](http://rp-bp.readthedocs.io/en/latest/) for instructions on how to install the package, or to use Docker/Singularity containers with the package pre-installed. Detailed usage instructions and tutorials are available.

## How to report issues

For bugs, issues, or feature requests, use the [bug tracker](https://github.com/dieterich-lab/rp-bp/issues). Follow the instructions and guidelines given in the templates.

## How to cite

Brandon Malone, Ilian Atanassov, Florian Aeschimann, Xinping Li, Helge Großhans, Christoph Dieterich. [Bayesian prediction of RNA translation from ribosome profiling](https://doi.org/10.1093/nar/gkw1350), _Nucleic Acids Research_, Volume 45, Issue 6, 7 April 2017, Pages 2960-2972.

## License

The MIT License (MIT). Copyright (c) 2016 dieterich-lab.
