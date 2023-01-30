# Ribosome profiling with Bayesian predictions (Rp-Bp)

**Rp-Bp** is an unsupervised Bayesian approach to predict translated open reading frames (ORFs) from ribosome profiles. **Rp-Bp** can be used for ORF discovery, or simply to estimate periodicity in a set of Ribo-seq samples.

**Rp-Bp** comes with two *interactive dashboards* or *web applications*, one for read and periodicity quality control, the other to facilitate Ribo-seq ORFs discovery. 

![rpbp](docs/images/logo-rpbp-final.png)

---

Ribosome profiling (Ribo-seq) is an RNA-sequencing-based readout of RNA translation. Isolation and deep-sequencing of ribosome-protected RNA fragments (ribosome footprints) provides a genome-wide snapshot of the translatome at sub-codon resolution. **Rp-Bp** is an unsupervised Bayesian approach to predict translated open reading frames (ORFs) from ribosome profiles, using the automatic Bayesian Periodic fragment length and ribosome P-site offset Selection (BPPS), *i.e.* read lengths and ribosome P-site offsets are inferred from the data, without supervision. **Rp-Bp** is able to handle *de novo* translatome annotation by directly assessing the periodicity of the Ribo-seq signal.

**Rp-Bp** can be used for ORF discovery, or simply to estimate periodicity in a set of Ribo-seq replicates, *e.g.* to know which samples and read lengths are usable for downstream analyses. When used for ORF discovery, **Rp-Bp** automatically classifies ORFs into different biotypes or categories, relative to their host transcript. 

## Documentation

Consult the [user guide](http://rp-bp.readthedocs.io/en/latest/) for instructions on how to install the package, or to use Docker/Singularity containers with the package pre-installed. Detailed usage instructions and tutorials are available. 

## How to report issues

Bugs and issues should be reported in the [bug tracker](https://github.com/dieterich-lab/rp-bp/issues). Follow the instructions and guidelines given in the template.

## How to cite

Brandon Malone, Ilian Atanassov, Florian Aeschimann, Xinping Li, Helge Großhans, Christoph Dieterich. [Bayesian prediction of RNA translation from ribosome profiling](https://doi.org/10.1093/nar/gkw1350), *Nucleic Acids Research*, Volume 45, Issue 6, 7 April 2017, Pages 2960-2972.

## License

The MIT License (MIT). Copyright (c) 2016 dieterich-lab.
