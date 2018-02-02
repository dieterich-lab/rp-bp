
# Running the Rp-Bp pipeline step-by-step

The Rp-Bp pipeline consists of an index creation step (refer to [Creating reference genome indices](#creating-reference-genome-indices)), which must be performed once for each genome and set of annotations, and a two-phase prediction pipeline, which must be performed for each sample (refer to [Running the Rp-Bp pipeline](#running-pipelines)). 

In the first phase of the prediction pipeline, a filtered genome profile is created (refer ro [Creating ORFs profiles](#creating-filtered-genome-profiles) for details). 

In the second phase, the ORFs which show evidence of translation in the profile are predicted (refer to [Predicting translated open reading frames](#predicting-translated-open-reading-frames) for details).

This document describes the steps to run each of these phases. It also shows some sample calls. For all programs, the `--help` option can be given to see the complete list of parameters.

<a id="toc"></a>

* [Creating reference genome indices](#creating-reference-genome-indices)
    * [Configuration file for index creation](#config-file-1)
    * [Input files](#creating-reference-genome-indices-input-files)
    * [Output files](#creating-reference-genome-indices-output-files)
    * [More about ORFs labels](#orfs-labels)
    * [More about *de novo* assembled transcripts](#de-novo-assembly-info)
* [Running the Rp-Bp pipeline](#running-pipelines)
    * [Configuration file for running the pipeline](#config-file-2)
    * [More about biological replicates](#more-about-replicates)
    * [1. Creating ORFs profiles](#creating-filtered-genome-profiles)
        * [1. More about the configuration file for running the pipeline](#config-file-more-1)
        * [1. Input files](#running-pipelines-input-1)
        * [1. Output files](#running-pipelines-output-1)
    * [2. Predicting translated ORFs](#predicting-translated-open-reading-frames)
        * [2. More about the configuration file for running the pipeline](#config-file-more-2)
        * [2. Input files](#running-pipelines-input-2)
        * [2. Output files](#running-pipelines-output-2)
    * [Using existing alignment files](#using-existing-alignment-files)
* [Logging options](#logging-options)
* [Parallel processing options](#parallel-processing-options) (SLURM options)

---

<a id="creating-reference-genome-indices"></a>

## Creating reference genome indices

This section describes the steps necessary to prepare a reference genome and
matching annotations for use in the Rp-Bp pipeline. The process must only be run
once for each reference genome and set of annotations.

The entire index creation process can be run automatically using the following command:

```
prepare-rpbp-genome <config> [--overwrite] [logging options] [processing options]
```

### Command line options

* `config` A [YAML](http://www.yaml.org/start.html) configuration file, as described below. A sample configuration file is also available with the example ([Download the dataset](#running-example)). 
* [`--overwrite`] Unless this flag is given, then steps for which the output files already exist will be skipped.
* [`logging options`] See [logging options](#logging-options).
* [`processing options`] See [parallel processing options](#parallel-processing-options).

<a id="config-file-1"></a>

### Configuration file for index creation

The following keys are read from the configuration file. Keys with [`brackets`] are optional. See also [Input files](#creating-reference-genome-indices-input-files) for format specifications.

* `gtf` The path to the reference annotations (format GTF2 or GFF3). 
* [`de_novo_gtf`] An additional GTF/GFF containing annotations constructed from a _de novo_ assembly (*e.g.* from StringTie). This file and the reference annotation file **must** have the same format (both GTF or GFF3). See also [More about *de novo* assembled transcripts](#de-novo-assembly-info).
* `fasta` The path to the reference genome sequence.
* `ribosomal_fasta` The path to the ribosomal sequence.

* `genome_base_path` The path to the output directory for the transcript fasta and ORFs.
* `genome_name` A descriptive name to use for the created files.
* `ribosomal_index` The base output path for the Bowtie2 index of the ribosomal sequence.
* `star_index` The base output path for the STAR index for the reference genome sequence. *The STAR index does not include the transcript annotations. They are incorporated on the fly during the mapping step.*

* [`orf_note`] An additional description used in the filename of the created ORFs.
* [`start_codons`] A list of strings that will be treated as start codons when searching for ORFs. Default: [`ATG`].
* [`stop_codons`] A list of strings that will be treated as stop codons when searching for ORFS. Default: [`TAA`, `TGA`, `TAG`].

<a id="creating-reference-genome-indices-input-files"></a>

### Input files

The required input files are those suggested by the configuration file keys described above. Please see [More about creating reference indices](#creating-reference-indices.md) for additional
details about choosing the appropriate files (*e.g.* from Ensembl).

* `gtf` The GTF/GFF file containing the reference annotations. The `exon` features must be present, and the transcript identifiers (`transcript_id` attribute for GTF) must match for exons from the same transcript. As the ORFs are labeled based on their positions relative to annotated coding sequences, the `CDS` features must also be present in the annotations. The `START codon` should be included in the CDS, but the `STOP codon` should not. *By default, Rp-Bp treats annotations according to the GTF2 (GTF) specifications, when the file extension is .gtf. Limited support is provided for [GFF3 format specifications](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md). When the file extension is .gff, it will be treated according to the GFF3 format specifications. In this case, it is assumed that the START and STOP codons are both included in the CDS. Rp-Bp will then automatically remove the STOP codon during processing.*
* [`de_novo_gtf`] The GTF/GFF file containing the _de novo_ assembled transcripts. It **must** have the same format specifications at the `gtf` file.

* `fasta` The input fasta file should contain all chromosome sequences. The identifiers must match those in the GTF/GFF file. Typically, the "primary assembly" file from Ensembl contains the appropriate sequences and identifiers. Please see [Ensembl](http://www.ensembl.org/info/genome/genebuild/assembly.html) for more information about the differences between assemblies.

* `ribosomal_fasta` The ribosomal DNA sequence is typically repeated many times throughout the genome. Consequently, it can be difficult to include it in the genome assembly and is often omitted. Therefore, a separate fasta file is required for these sequences (which are later used to filter reads). This file could also include other sequences which should be filtered, such as *e.g.* tRNA sequences downloaded from [GtRNAdb](http://gtrnadb.ucsc.edu/).

<a id="creating-reference-genome-indices-output-files"></a>

### Output files

The base path for these files is: `<genome_base_path>/`:

* `<ribosomal_index>/` The Bowtie2 index files (`<ribosomal_index>`.1.bt2, *etc*.) for the ribosomal_fasta file.

* `<star_index>/` The STAR index files (SA, Genome, *etc*.) for the `fasta` file.

* `<genome_name>.annotated.bed.gz` A BED file containing all transcripts, with the following columns:
    * *seqname*: name of the chromosome or scaffold.
    * *start, end*: start, end position of the transcript in standard chromosomal coordinates.
    * *id*: transcript id.
    * *score*: currently unused.
    * *strand*: defined as + (forward) or - (reverse).
    * *thick_start, thick_end*: coordinate at which to start drawing the CDS (-1 for noncoding transcripts).
    * *color*: currently unused.
    * *num_exons*: number of exons pertaining to the transcript.
    * *exon_lengths*: size of the exons.
    * *exon_genomic_relative_starts*: start coordinate of each exon.

#### Output from the annotations

The base path for these files is: `<genome_base_path>/transcript-index/`:

   * `<genome_name>.transcripts.annotated.fa`. The sequences of all annotated transcripts.
    
   * `<genome_name>.genomic-orfs.annotated.<orf_note>.bed.gz`. A BED file containing all ORFs from annotated transcripts. In addition to the standard BED columns as described above (but with respect to ORFs), this file also includes columns giving the `orf_type`, `orf_length`, and `orf_num`. The ORF `id`s are of the form: `transcript_seqname:start-end:strand`. The START codon _is_ included in the ORF, but the STOP codon _is not_. The `thick_start` and `thick_end` are always the same as `start` and `end`.
    
   * `<genome_name>.orfs.annotated.<orf_note>.bed.gz`. A BED file containing each exon from each ORF. The `id` of an exon corresponds exactly to the ORF to which it belongs; exons from the same ORF have the same `id`. The extra columns are `exon_index`, which gives the order of the exon in the transcript, and `transcript_start`, which gives the start position of that index in transcript coordinates. Exons are always sorted by "lowest start first", so the order is really reversed for ORFs on the reverse strand.
  
#### Output from a _de novo_ assembly 

See also [More about *de novo* assembled transcripts](#de-novo-assembly-info). The semantics of these files are the same as above, but they are created using the `de_novo_gtf` annotations. 
*N.B. ORFs which completely overlap annotations are not included.*

The base path for these files is: `<genome_base_path>/transcript-index/`:

  * `<genome_name>.transcripts.de-novo.fa`
  * `<genome_name>.genomic-orfs.de-novo.<orf_note>.bed.gz`
  * `<genome_name>.orfs.de-novo.<orf_note>.bed.gz`

#### Output from both the annotations and a _de novo_ assembly 

The semantics are again the same as above. If a _de novo_ assembly was not provided, these are simply symlinks to the respective "annotations" files. Otherwise, they are the concatenation of the respective "annotation" and "_de novo_" files.

  * `<genome_base_path>/<genome_name>.star-input.gtf`
  
The base path for these files is: `<genome_base_path>/transcript-index/`:
 
  * `<genome_name>.genomic-orfs.<orf_note>.bed.gz`
  * `<genome_name>.orfs.<orf_note>.bed.gz`

<a id="orfs-labels"></a>

### More about ORFs labels

As part of the index creation phase, ORFs are labeled according to their location relative to the annotated, canonical transcripts, *i.e.* the annotated CDS regions. We use the following labels (essentially, the same as those given in the supplement of the paper):

* `canonical`: an ORF which exactly coincides with an annotated transcript.
* `canonical_extended`: an ORF which starts upstream of a `canonical` ORF, but otherwise has the same structure.
* `canonical_truncated`: an ORF which starts downstream of a `canonical` ORF, but otherwise has the same structure.
* `five_prime`: an ORF which is completely in the annotated 5' leader region of a transcript and does not overlap other `canonical` ORFs on the same strand (such as from another isoform).
* `three_prime`: an ORF in the annotated 3' trailer region of a transcript and does not overlap other `canonical` ORFs on the same strand.
* `noncoding`: an ORF from a transcript not annotated as coding, such as a lncRNA or pseudogene.
* `five_prime_overlap`: an ORF in the annotated 5' leader region of a transcript but which overlaps a `canonical` ORF; this could either be out-of-frame with respect to the `canonical` ORF of the same transcript or otherwise overlap a `canonical` ORF from another isoform.
* `three_prime_overlap`: an ORF in the annotated 3' trailer region of a transcript but which overlaps a `canonical` ORF.
* `suspect`: an ORF which partially overlaps the interior of a `canonical` ORF. These can result from things like retained introns.
* `within`: an out-of-frame ORF in the interior of a `canonical` ORF.

<a id="de-novo-assembly-info"></a>

### More about _de novo_ assembled transcripts

In principle, there is no difference between "standard" annotated transcripts and _de novo_ assembled transcripts. In both cases, ORFs are extracted from the transcripts based on the given start and stop codons. However, it is often of scientific interest to distinguish between ORFs from annotated transcripts and "novel" ones from a _de novo_ assembly.

Both to avoid reporting spurious "novel" translated ORFs and to avoid redundant calculations, we remove each ORF identified in the _de novo_ assembly which exactly match an ORF using only the annotations.

We then label ORFs from the _de novo_ assembly as follows, with respect to their relationship to the annotated transcripts. Due to the first filtering step, the following categories include only ORFs which are at least partially only supported in the _de novo_ assembly.

* `novel`: The ORF does not overlap the annotations at all.
* `novel_canonical_extended`: an ORF which starts upstream of a `canonical` ORF, but otherwise has the same structure.
* `novel_noncoding`: an ORF from a transcript not annotated as coding, such as a lncRNA or pseudogene.
* `novel_five_prime_overlap`: an ORF in the annotated 5' leader region of a transcript but which overlaps a `canonical` ORF; this could either be out-of-frame with respect to the `canonical` ORF of the same transcript or otherwise overlap a `canonical` ORF from another isoform.
* `novel_three_prime_overlap`: an ORF in the annotated 3' trailer region of a transcript but which overlaps a `canonical` ORF.
* `novel_suspect`: an ORF which partially overlaps the interior of a `canonical` ORF. These can result from things like retained introns.
* `novel_within`: an out-of-frame ORF in the interior of a `canonical` ORF

The other labels, such as `novel_canonical` or `novel_five_prime` are guaranteed to be filtered in the first step based on their definitions.

In preparation to running the main pipeline, `reference` and `de_novo` annotations are concatenated in a single file which is used on the fly by `STAR` at the mapping step. Merging and sorting the resulting annotations is not necessary.

[Back to top](#toc)

<a id='running-pipelines'></a>

## Running the Rp-Bp pipeline

In the first phase of the prediction pipeline, a filtered genome profile is created, as explained in [Creating ORF profiles](#creating-filtered-genome-profiles).

In the second phase, the ORFs which show evidence of translation in the profile are predicted, as explained in [Predicting translated ORFs](#predicting-translated-open-reading-frames).

[More details about biological replicates](#more-about-replicates) are given below.


The entire Rp-Bp pipeline can be run on a set of riboseq samples (including any biological replicates) which all use the same genome indices using a sample sheet-like [YAML](http://www.yaml.org/start.html) configuration file with the following command:

```
# Do not merge the replicates.
run-all-rpbp-instances <config> [--tmp] [--overwrite] [-k/--keep-intermediate-files] [--flexbar-options] [--star-executable] [--star-read-files-command] [--star-additional-options] [logging options] [processing options]

# Merge the replicates, do not calculate Bayes factors nor make predictions for individual samples.
run-all-rpbp-instances <config> --merge-replicates [--tmp] [--overwrite] [-k/--keep-intermediate-files] [--flexbar-options] [--star-executable] [--star-read-files-command] [--star-additional-options] [logging options] [processing options]

# Merge the replicates and also calculate Bayes factors and make predictions for individual samples.
run-all-rpbp-instances <config> --merge-replicates --run-replicates [--tmp] [--overwrite] [-k/--keep-intermediate-files] [--flexbar-options] [--star-executable] [--star-read-files-command] [--star-additional-options] [logging options] [processing options]
``` 
    
### Command line options

* `config` The [YAML](http://www.yaml.org/start.html) configuration file, as described below. A sample configuration file is also available with the example ([Download the dataset](#running-example)). 
* [`--tmp <loc>`] If this flag is given, then all relevant calls will use `<loc>` as the base temporary directory. Otherwise, the program defaults will be used.
* [`--overwrite`] Unless this flag is given, then steps for which the output files already exist will be skipped.
* [`-k/--keep-intermediate-files`] Unless this flag is given, large intermediate files, such as fastq files output by flexbar after removing adapters, will be deleted.
* [`--flexbar-options`] A space-delimited list of options to pass to flexbar. Each option must be quoted separately as in "--flexbarOption value". For quality-based trimming *e.g.*, one may pass the quality-based trimming mode and format. Default: see [Creating filtered genome profiles](#creating-filtered-genome-profiles).
* [`--star-executable`] In principle, `STARlong` (as opposed to `STAR`) could be used for alignment. Given the nature of riboseq reads (that is, short due to the experimental protocols of degrading everything not protected by a ribosome), this is unlikely to be a good choice, though. Default: `STAR`.
* [`--star-read-files-command`] The input for `STAR` will always be a gzipped fastq file. `STAR` needs the system command which means "read a gzipped text file". The program attempts to guess the name of this command based on the operating system (*e.g.* OSX, Ubuntu), but it can be explicitly specified as a command line option. Default: `gzcat` if `sys.platform.startswith("darwin")`; `zcat` otherwise. Please see [python.sys documentation](https://docs.python.org/3/library/sys.html) for more details about attempting to guess the operating system.
* [`--star-additional-options`] A space-delimited list of additional options to pass to star. Each option must be quoted separately as in "--starOption value". Default: see [Creating filtered genome profiles](#creating-filtered-genome-profiles).
* [`logging options`] See [logging options](#logging-options).
* [`processing options`] See [parallel processing options](#parallel-processing-options).

<a id="config-file-2"></a>

### Configuration file for running the pipeline

The following keys are read from the configuration file. Keys with [`brackets`] are optional.
Additional details about the configuration file options are given in *More about the configuration file for running the pipeline* [(part 1)](config-file-more-1) and [(part 2)](config-file-more-2).

* `riboseq_data` The base output location for all created files.
* [`note`] An optional string which will be added to all file names. It should not contain spaces or any other special characters.
* [`models_base`] The base path to the compiled models. The models specified in the paper are included with the source distribution and compiled/pickled as part of the installation process. They are installed in an operating system-specific location (in particular, `user_data_dir` from the [appdirs package](https://pypi.python.org/pypi/appdirs). This location is determined during installation and does not normally need to be changed. For development, they may be in some alternative location.

* `riboseq_samples` A dictionary in which each entry specifies a sample. The key is an informative name about the sample, and the value gives the complete path to the sequencing file (a fastq or fastq.gz file). The names will be used to construct filenames, so they should not contain spaces or other special characters.
* `riboseq_biological_replicates` A dictionary in which each entry species one condition and all samples which are replicates of the condition. The key of the dictionary is a string description of the condition, and the value is a list that gives all of the sample replicates which belong to that condition. The names of the sample replicates must match the dataset names specified in `riboseq_samples`.

In addition, the options below should be exactly the same as those used in the [Configuration file for index creation](#config-file-1). In practice, the same configuration file can be used for index creation and running the pipeline.

* `gtf` The path to the reference annotations.
* `fasta` The path to the reference genome sequence.
* `genome_base_path` The path to the output directory for the transcript fasta and ORFs.
* `genome_name` A descriptive name to use for the created files.
* `ribosomal_index` The base path for the Bowtie2 index of the ribosomal sequence.
* `star_index` The base path to the STAR index.

<a id="more-about-replicates"></a>

### More about biological replicates

The Rp-Bp pipeline handles replicates by adding the (smoothed) ORF profiles. The Bayes factors and predictions are then calculated based on the combined profiles. The `--merge-replicates` flag indicates that the replicates should be merged. By default, if the `--merge-replicates` flag is given, then predictions will not be made for the individual datasets, unless the `--run-replicates` flag is also given, in which case predictions will be made for both the merged replicates as well as the individual samples.

[Back to top](#toc)

<a id='creating-filtered-genome-profiles'></a>

### 1. Creating ORF profiles

The entire profile creation process can be run automatically using the `create-orf-profiles` script. This script is called by default when running the main pipeline (when calling `run-all-rpbp-instances`). If one is interested in only creating the filtered genome profiles, then `create-orf-profiles` can be called separately:

```
create-orf-profiles <raw_data> <config> <sample name> [--tmp] [--overwrite] [-k/--keep-intermediate-files] [--flexbar-options] [--star-executable] [--star-read-files-command] [--star-additional-options] [logging options] [processing options]
```

#### Command line options

* `config` The [YAML](http://www.yaml.org/start.html) configuration file, as described below. *The script reads most of the required paths from the configuration file, so if running separately, the arguments must be consistent with the paths given in the configuration file.*
* [`--tmp <loc>`] If this flag is given, then all relevant calls will use `<loc>` as the base temporary directory. Otherwise, the program defaults will be used.
* [`--overwrite`] Unless this flag is given, then steps for which the output files already exist will be skipped.
* [`-k/--keep-intermediate-files`] Unless this flag is given, large intermediate files, such as fastq files output by flexbar after removing adapters, will be deleted.
* [`--flexbar-options`] A space-delimited list of options to pass to flexbar. Each option must be quoted separately as in `--flexbarOption value`. For quality-based trimming *e.g.*, one may pass the quality-based trimming mode and format. Default: see below.
* [`--star-executable`] In principle, `STARlong` (as opposed to `STAR`) could be used for alignment. Given the nature of riboseq reads (that is, short due to the experimental protocols of degrading everything not protected by a ribosome), this is unlikely to be a good choice, though. Default: `STAR`.
* [`--star-read-files-command`] The input for `STAR` will always be a gzipped fastq file. `STAR` needs the system command which means "read a gzipped text file". The program attempts to guess the name of this command based on the operating system (*e.g.* OSX, Ubuntu), but it can be explicitly specified as a command line option. Default: `gzcat` if `sys.platform.startswith("darwin")`; `zcat` otherwise. Please see [python.sys documentation](https://docs.python.org/3/library/sys.html) for more details about attempting to guess the operating system.
* [`--star-additional-options`] A space-delimited list of additional options to pass to star. Each option must be quoted separately as in `--starOption value`. Default: see below.
* [`logging options`] See [logging options](#logging-options).
* [`processing options`] See [parallel processing options](#parallel-processing-options).

<a id='config-file-more-1'></a>

#### 1. More about the configuration file for running the pipeline

The following same keys are read from the configuration file, see [Configuration file for running the pipeline](#config-file-2).

* `riboseq_data` The base output location for all created files.
* [`note`] An optional string which will be added to all file names. It should not contain spaces or any other special characters.
* [`models_base`] The base path to the compiled models.

* `gtf` The path to the reference annotations.
* `genome_base_path` The path to the output directory for the transcript fasta and ORFs.
* `genome_name` A descriptive name to use for the created files.
* `ribosomal_index` The base path for the Bowtie2 index of the ribosomal sequence.
* `star_index` The base path to the STAR index.

##### Flexbar and STAR options

A number of pre-defined options can be passed to flexbar and STAR via the configuration file.
Any additional options need to be passed to flexbar and STAR as arguments in the same way as when calling `run-all-rpbp-instances`, see [Running the Rp-Bp pipeline](running-pipelines).
If the same option is used in both locations (in the configuration file and passed as an argument, then the value passed as an argument has precedence).

**N.B.** Keys in the configuration file do NOT correspond to the name of the respective options used by flexbar or STAR. However, when passing options as arguments, one must use the exact semantic required by the external program in question (there is no "validity check"). Note also that some options for STAR are handled using `--num-cpus` or `--mem`, which cannot be overridden, see [Parallel processing options](#parallel-processing-options).

**N.B.** These do not normally need to be given, unless the default values are to be modified.

###### Flexbar options
* [`adapter_file`] A fasta file containing a set of adapter sequences used by flexbar. Recommend via configuration file. Default: None.
* [`adapter_sequence`] A single sequence used to remove adapters within flexbar. Default: None.
* [`max_uncalled`] The maximum number of Ns to permit in a read without filtering. Default 1.
* [`pre_trim_left`] The number of bases to remove from the 5' end of all reads: Default: 0.

###### STAR options
* [`align_intron_min`] Default: 20.
* [`align_intron_max`] Default: 100000.
* [`out_filter_intron_motifs`] Default: RemoveNoncanonicalUnannotated.
* [`out_filter_mismatch_n_max`] Default: 1.
* [`out_filter_mismatch_n_over_l_max`] Default: 0.04.
* [`out_filter_type`] Default: BySJout.
* [`out_sam_attributes`] Default: AS NH HI nM MD.
* [`sjdb_overhang`] Default: 50.

##### Rp-Bp specific options (via the configuration file only)

**N.B.** These do not normally need to be given, unless the default values are to be modified.

###### Multimapper options
* [`keep_riboseq_multimappers`] If this variable is present in the configuration file with any value (even something like "no" or "null" or "false"), then multimapping riboseq reads *will not* be removed. They be treated as "normal" reads in every place they map, *i.e.* the weight of the read will not be distributed fractionally, probabilistically, *etc.*

###### Metagene periodicity options
* [`seqids_to_keep`] If this list is given, then only transcripts appearing on these identifiers will be used to construct the metagene profiles (and other downstream analysis). The identifiers must match exactly (*e.g.* "2" and "chr2" do not match).

* [`metagene_profile_start_upstream`] The number of bases upstream of the translation initiation site to begin constructing the metagene profile. Default: 50.
* [`metagene_profile_start_downstream`] The number of bases downstream of the translation initiation site to end the metagene profile. Default: 20.
* [`metagene_profile_end_upstream`] The number of bases upstream of the translation termination site to begin constructing the metagene profile. Default: 50.
* [`metagene_profile_end_downstream`] The number of bases downstream of the translation termination site to end the metagene profile. Default: 20.

* [`periodic_offset_start`] The position, relative to the translation initiation site, to begin calculating periodicity Bayes factors (inclusive)
* [`periodic_offset_end`] The position, relative to the translation initiation site, to stop calculating periodicity Bayes factors (inclusive)
* [`metagene_profile_length`] The length of the profile to use in the models. `metagene_profile_length` + `periodic_offset_end` must be consistent with the length of the extracted metagene profile

* [`metagene_profile_iterations`] The number of iterations to use for each chain in the MCMC sampling. The first half of the iterations are discarded as burn-in samples. All of the remaining samples are used to estimate the posterior distributions. That is, we do not use thinning. default: 500

###### Periodicity and offset options
* [`use_fixed_lengths`] If this variable is present in the config file with any value (even something like "no" or "null" or "false"), then the estimated  periodic read lengths and offsets will not be used. Instead, fixed values given by `lengths` and `offsets` (below) will be used.
* [`lengths`] A list of read lengths which will be used for creating the profiles if the `use_fixed_lengths` option is given. Presumably, these are lengths that have periodic metagene profiles.
* [`offsets`] The P-site offset to use for each read length specifed by `lengths` if the `use_fixed_lengths` option is given. The number of offsets must match the number of lengths, and they are assumed to match. For example `lengths` of 26, 29 with `offsets` 9, 12 means only reads of lengths 26 bp and 29 bp will be used to create the profiles. The 26 bp reads will be shifted by 9 bp in the 5' direction, while reads of length 29 bp will be shifted by 12 bp.

* [`min_metagene_profile_count`] If fixed lengths are NOT used: the minimum number of reads for a particular length in the filtered genome profile. Read lengths with fewer than this number of reads will not be used. Default: 1000.
* [`min_metagene_bf_mean`] If fixed lengths are NOT used: if `max_metagene_profile_bayes_factor_var` is not `None`, then this is taken as a hard threshold on the estimated Bayes factor mean. Default: 5.
* [`max_metagene_bf_var`] If fixed lengths are NOT used: if given, then this is taken as a hard threshold on the estimated Bayes factor variance. Default: None, *i.e.* this filter is not used. (`null` in yaml)
* [`min_metagene_bf_likelihood`] If fixed lengths are NOT used: if given, then this is taken a threshold on the likelihood of periodicity. Default: 0.5.

If `min_metagene_bf_likelihood` is given, then this is taken as the boundary value; *i.e* a profile is "periodic" if:

```
[P(bf > min_metagene_bf_mean)] > min_metagene_bf_likelihood
```

If both `max_metagene_bf_var` and `min_metagene_bf_likelihood` are `None`, then this is taken as a hard threshold on the mean for selecting periodic read lengths.
If both `max_metagene_bf_var` and `min_metagene_bf_likelihood` are given, then both filters will be applied and the result will be the intersection. 

###### Smoothing options

* [`smoothing_fraction`] The fraction of the profile to use for smoothing within LOWESS. Default: 0.2.
* [`smoothing_reweighting_iterations`] The number of reweighting iterations to use within LOWESS. Please see the statsmodels documentation for a detailed description of this parameter. Default: 0.

* [`min_orf_length`] If this value is greater than 0, then ORFs whose length (in nucleotides) is less than this value will not be smoothed, and neither the Bayes factor estimates (nor the chi-square p-value) will be calculated. Default: 0.
* [`max_orf_length`] If this value is greater than 0, then ORFs whose length (in nucleotides) is greater than this value will not be smoothed, and neither the Bayes factor estimates (nor the chi-square p-value) will be calculated. Default: 0.
* [`min_signal`] ORFs for which the number of in-frame reads is less than this value will not be smoothed, and neither the Bayes factor estimates nor the chi-square p-value will be calculated. Default: 5.

###### Shared MCMC options
These affect the MCMC both for estimating metagene profile periodicity and ORF translation Bayes factors.

* [`seed`] The random seed for the MCMC sampling. Default: 8675309.
* [`chains`] The number of chains to use in the MCMC sampling. Default: 2.

<a id='running-pipelines-input-1'></a>

#### 1. Input files

* All the input files are those suggested by the configuration file keys, as explained above. 
* The models (`.pkl` files), located in `<models_base>/periodic` and `<models_base>/nonperiodic`.
* In addition, the input annotations `<genome_base_path>/<genome_name>.star-input.gtf` are used by `STAR`. For details, see [Creating reference genome indices](#creating-reference-genome-indices).

<a id='running-pipelines-output-1'></a>

#### 1. Output files

**N.B.** If the `keep_riboseq_multimappers` configuration option is given, then the `-unique` part will not be present in the output file names.

* Base path for trimmed and quality filtered reads: `<riboseq_data>/without-adapters/`
    * **trimmed and filtered reads** A fastq.gz file containing the reads after removing adapters and low-quality reads. 
    `<sample-name>[.<note>].fastq.gz`    
* Base path for reads aligning to ribosomal sequences: `<riboseq_data>/with-rrna/`
    * **discarded reads** A fastq.gz file containing reads which align to the ribosomal index. They are recorded but not used in later processing. 
    `<sample-name>[.<note>].fastq.gz`
* Base path for reads not aligning to ribosomal sequences: `<riboseq_data>/without-rrna/`
    * **retained reads** A fastq.gz file containing reads which do not align to the ribosomal index and are used in further processing. 
    `<sample-name>[.<note>].fastq.gz`
* Base path for aligned reads: `<riboseq_data>/without-rrna-mapping/`
    * **sorted reads aligned to the genome** A sorted bam file containing all alignments of reads to the genome. 
    `<sample-name>[.<note>]Aligned.sortedByCoord.out.bam`. 
    `<sample-name>[.<note>].bam`, which is a symlink to the `Aligned.sortedByCoord.out.bam` file.
    * **aligned reads which map uniquely to the genome** A sorted bam file containing all alignments of reads to the genome with multimapping reads filtered out. 
    `<sample-name>[.<note>]-unique.bam`

Indices are also created for the bam files. In addition, `STAR` creates some files in the location `<riboseq_data>/without-rrna-mapping/<sample-name>_STARgenome`.

* Base path for metagene profiles: `<riboseq_data>/metagene-profiles/`
    * **metagene profiles** A gzipped csv file containing the metagene profiles (given by the `position` or offset and `count` columns) for all read lengths (`length` column) which occur in the uniquely-aligning reads. It includes the metagene profile around both the annotated translation initiation site and translation termination site (`type` column). 
    `<sample-name>[.<note>]-unique.metagene-profile.csv.gz`
    * **periodicity estimations** A gzipped csv file containing the Bayes factor estimates for all P-site offsets found or specified, as well as information about the number of reads in the respective profile. 
    `<sample-name>[.<note>]-unique.metagene-periodicity-bayes-factors.csv.gz`
    * **estimated P-site offsets** A gzipped csv file containing the selected P-site offset for each read length. All read lengths are included, even if the estimates do not meet the criteria specified in the configuration file. (The filtering occurs later.) 
    `<sample-name>[.<note>]-unique.periodic-offsets.csv.gz`
* Base path for ORF profiles: `<riboseq_data>/orf-profiles/`
    * **unsmoothed ORF profiles** A gzipped, sparse [matrix market file](http://math.nist.gov/MatrixMarket/formats.html) containing the profiles for all ORFs (`orf_num`, `orf_position` and `read_count`). **N.B.** The matrix market format uses base-1 indices. 
    `<sample-name>[.<note>]-unique.length-<lengths>.offset-<offsets>.profiles.mtx.gz` 

The smoothed profiles are not explicitly stored. 

#### Difference from paper

The fifth step of creating the base genome profile in the paper is "Everything except the 5' end of the remaining reads is removed." This profile is not explicitly constructed in the pipeline. The `<sample_name>-unique.bam` file already contains the necessary information.

[Back to top](#toc)

<a id='predicting-translated-open-reading-frames'></a>

### 2. Predicting translated ORFs

The entire translation prediction process can be run automatically using the `predict-translated-orfs` script. This script is called by default when running the main pipeline (when calling `run-all-rpbp-instances`). If one is interested in translation prediction, given that ORFs profiles are already available, then `predict-translated-orfs` can be called separately:

```
predict-translated-orfs <config> <sample or condition name> [--overwrite] [--merge-replicates] [logging options] [processing options]
```

#### Command line options

* `config` The [YAML](http://www.yaml.org/start.html) configuration file, as described below. *The script reads most of the required paths from the configuration file, so if running separately, the arguments must be consistent with the paths given in the configuration file.*
* `sample or condition name` The name of either one of the `riboseq_samples` or `riboseq_biological_replicates` from the configuration file (if merging replicates).
* [`--overwrite`] Unless this flag is given, then steps for which the output files already exist will be skipped.
* [`--merge-replicates`] If this flag is present, then the ORF profiles will be merged for all replicates in the condition given by `<sample or condition name>`. If this flag is is present, the `--overwrite` flag will automatically be set.
* [`logging options`] See [logging options](#logging-options).
* [`processing options`] See [parallel processing options](#parallel-processing-options).

<a id='config-file-more-2'></a>

#### 2. More about the configuration file for running the pipeline

The following same keys are read from the configuration file, see [Configuration file for running the pipeline](#config-file-2).

* `riboseq_data` The base output location for all created files.
* [`note`] An optional string which will be added to all filen ames. It should not contain spaces or any other special characters.
* [`models_base`] The base path to the compiled models.

* `fasta` The path to the reference genome sequence.
* `genome_base_path` The path to the output directory for the transcript fasta and ORFs.
* `genome_name` A descriptive name to use for the created files.

##### Rp-Bp specific options (via the configuration file only)

**N.B.** These do not normally need to be given, unless the default values are to be modified.

###### Bayes factor estimate options
* [`chi_square_only`] If this flag is in the config file with any value, then only the Rp-chi pipeline will be performed; namely, the translation models will not be fit to the data, and the posterior distributions will not be estimated. Otherwise, only th e Rp-Bp pipeline is run.

* [`translation_iterations`] The number of iterations to use for each chain in the MCMC sampling. The first half of the iterations are discarded as burn-in samples. All of the remaining samples are used to estimate the posterior distributions. That is, we do not use thinning. Default: 200.

###### Selecting predicted ORFs options
* [`min_bf_mean`] The minimum value for the estimated Bayes factor mean to "predict" that an ORF is translated. This value is used in conjunction with both `min_bf_mean` and `min_bf_likelihood`. Default: 5.
* [`max_bf_var`] The maximum value value for the estimated Bayes factor variance to "predict" that an ORF is translated. ORFs must meet both the `min_bf_mean` and `max_bf_var` filters to be predicted. If `max_bf_var` is a positive value, then this is taken as a hard threshold on the estimated Bayes factor mean. ORFs must meet both the `min_bf_mean` and `max_bf_var` filters to be selected as "translated." Default: null (*i.e.* this filter is not used by default).
* [`min_bf_likelihood`] The minimum probability of the BF exceeding `min_bf_mean` to select an ORF as translated. Default: 0.5.

If `min_bf_likelihood` is given, then this is taken as the boundary value; *i.e.* an ORF is selected as "translated" if:

```
[P(bf > min_bf_mean)] > min_bf_likelihood
```

If both `max_bf_var` and `min_bf_likelihood` are `None` (null in YAML), then this is taken as a hard threshold on the mean for selecting translated ORFs.
If both `max_bf_var` and `min_bf_likelihood` are given, then both filters will be applied and the result will be the intersection.
 
* [`chisq_significance_level`] Only used with `chi_square_only`. For the chi-square test, this value is first Bonferroni corrected based on the number of ORFs which pass the smoothing filters. It is then used as the significance threshold to select translated ORFs. Default: 0.01.

###### Shared MCMC options
These affect the MCMC both for estimating metagene profile periodicity and ORF translation Bayes factors, see [1. More about the configuration file for running the pipeline](#config-file-more-1).

* [`seed`] The random seed for the MCMC sampling. Default: 8675309.
* [`chains`] The number of chains to use in the MCMC sampling. Default: 2.

<a id='running-pipelines-input-2'></a>

#### 2. Input files

This script requires several files created during the previous steps of the pipeline, as well as a few external files. These would normally be given by the configuration file keys, as explained above, and are thus readily available when running the main pipeline (`run-all-rpbp-instances`).

* External files:
    * **genome fasta file** The genome fasta file. This is the same file used for `prepare-rpbp-genome`.
    * **orfs** The ORFs (gzipped BED file) created by `prepare-rpbp-genome`. 
    It must be located at `<genome_base_path>/transcript-index/<genome_name>.genomic-orfs.<orf_note>.bed.gz`
    * **exons** The ORF exons (gzipped BED file) created by `prepare-rpbp-genome`. 
    It must be located at `<genome_base_path>/transcript-index/<genome_name>.orfs.<orf_note>.bed.gz`
    * **models of translation** Some compiled, pickled Stan model files must be located in the `<models_base>/translated` folder.
    * **models of lack of translation** Some compiled, pickled Stan model files must be located in the `<models_base>/untranslated` folder.

* Metagene and ORF profiles: These files are the output files from the preceding phase, see [1. Output files](#running-pipelines-output-1). 
    
    The base directory for the metagene profiles-related files is: `<riboseq_data>/metagene-profiles/`:
    * **metagene profiles** `<sample-name>[.<note>]-unique.metagene-profile.csv.gz`    
    * **periodicity estimations** `<sample-name>[.<note>]-unique.metagene-periodicity-bayes-factors.csv.gz`
    * **estimated P-site offsets** `<sample-name>[.<note>]-unique.periodic-offsets.csv.gz`
    
    The base directory for the ORF profiles is: `<riboseq_data>/orf-profiles/`  
    * **unsmoothed ORF profiles** `<sample-name>[.<note>]-unique.length-<lengths>.offset-<offsets>.profiles.mtx.gz`
    
<a id='running-pipelines-output-2'></a>

#### 2. Output files
    
If replicates are merged, then these files will be created for each condition. Otherwise, they will be created for each sample (or both if the appropriate options are given). 
Furthermore, there are `unfiltered` and `filtered` versions of the files. The `filtered` version result from performing the filtering described in the paper (taking the longest predicted ORF for each stop codon, and then selecting the ORF with the highest expected Bayes factor among each group of overlapping ORFs). The `unfiltered` version contains all predictions.

The base path for all output files is: `<riboseq_data>/orf-predictions/`:

* Estimates
    * **the Bayes factor estimates** A BED file which contains the estimated values for all ORFs. The first 12 columns are valid BED12 entries that are simply copied from the `orfs` BED file (see [Output files](#creating-reference-genome-indices-output-files)). 
    `<sample_name>[.<note>]-unique.length-<lengths>.offset-<offsets>.bayes-factors.bed.gz` 
    
* Predictions
    * **Rp-Bp predictions** are made using the methodology described in the paper.
        * **the ORFs** A BED file containing the ORFs in the **final prediction** set. The format is identical to the "Bayes factor estimates" file above. 
        `<sample_name>[.<note>]-unique.length-<lengths>.offset-<offsets>.predicted-orfs.bed.gz`
        * **the DNA sequences** A fasta file containing the DNA sequences of the predicted ORFs. The fasta header matches the `id` column in the BED files. 
        `<sample_name>[.<note>]-unique.length-<lengths>.offset-<offsets>.predicted-orfs.dna.gz`    
        * **the protein sequences** A fasta file containing the protein sequences of the predicted ORFs. The fasta header matches the `id` column in the BED files. 
        `<sample_name>[.<note>]-unique.length-<lengths>.offset-<offsets>.predicted-orfs.protein.gz`
 
N.B. If `smoothing options` are specified in the configuration file with any value, the following string `<frac-fraction>.<rw-reweighting_iterations>` is also appended to the file names. For default values (unless they are explicitly given in the configuration file), this is omitted. 
  
If using the Rp-chi pipeline, then predictions are made using a simple chi-square test. The semantics of these files are the same as above.

N.B. These files will only be generated with the `chi_square_only` option.

   * **the ORFs**. `<sample_name>[.<note>]-unique.length-<lengths>.offset-<offsets>.chisq.predicted-orfs.bed.gz`
   * **the DNA sequences**. `<sample_name>[.<note>]-unique.length-<lengths>.offset-<offsets>.chisq.predicted-orfs.dna.gz`
   * **the protein sequences**. `<sample_name>[.<note>]-unique.length-<lengths>.offset-<offsets>.chisq.predicted-orfs.protein.gz`

[Back to top](#toc)

<a id='using-existing-alignment-files'></a>

### Using existing alignment files

While we found the `flexbar`-`bowtie2`-`STAR` pipeline to be effective for
processing and aligning reads, the probabilistic models of Rp-Bp do not depend
on any particular properties of these programs. Thus, it is possible to use any
alignment files for the pipeline. Essentially, the alignment (bam) files need
to be named as expected by the pipeline and placed in the correct folders.

Please see the [custom alignment files](#custom-alignment-files.md) documentation
for more details.

[Back to top](#toc)

<a id='logging-options'></a>

## Logging options

All of the driver scripts mentioned above, and many of the internal scripts as well, allow detailed specification of logging options on the command line. Interally, logging is handled using the standard python logging system. The following options are allowed.

The allowed logging levels are: `NOTSET`, `DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`. Most of the scripts output possibly-useful information at the `INFO` level; some scripts also output `DEBUG` messages, but as the name indicates, they are typically only useful for tracking down problems which do not cause actual errors.

* `--log-file` This option specifies a file to which logging statements will be written (in addition to stdout and  stderr, if specified). Default: False

* `--log-stdout` If this flag is present, then logging statements will be written to stdout (in addition to a file and stderr, if specified). Default: False

* `--log-stderr` If this flag is present, then logging statements will be written to stderr (in addition to a file and stdout, if specified). Default: True

* `--logging-level` If this value is specified, then it will be used for all logs. default: WARNING

* `--file-logging-level` The logging level to be used for the log file, if specified. This option overrides `--logging-level`. Default: NOTSET

* `--stdout-logging-level` The logging level to be used for the stdout log, if specified. This option overrides `--logging-level`. Default: NOTSET

* `--stderr-logging-level` The logging level to be used for the stderr log, if specified. This option overrides `--logging-level`. Default: NOTSET

[Back to top](#toc)

<a id='parallel-processing-options'></a>

## Parallel processing options

All of the scripts accept options for running code in parallel. Furthermore, the scripts are designed to seamlessly integrate with the [SLURM](http://slurm.schedmd.com/) scheduler. The following options are allowed. SLURM-specific options are specified by [`brackets`].

* `--num-cpus` The number of CPUs to use. The definition of a "CPU" varies somewhat among the programs. For example, for STAR, these are actually threads. For many of the python scripts, this number is translated into the number of processes to spawn. None of the code parallelizes across machines, so the value should not be greater than the number of cores on the machine on which the programs are executed. When used with SLURM, this will be translated into an sbatch request like: `--ntasks 1 --cpus-per-task <num-cpus>`. Default: 1

* `--mem` For STAR genome indexing, the amount of RAM to request. The rest of the programs do not use this value. When used with SLURM, this will be translated into an sbatch request like: `--mem=<mem>`. Default: 10G

* `--do-not-call` If this flag is present, then the commands will not be executed (but will be printed). This can be useful to ensure paths in configuration files are correct.

* [`--use-slurm`] If this flag is present, then the commands will be submitted to SLURM via sbatch. default: By default, each command is executed sequentially within the current terminal.

* [`--time`] The amount of time to request. This will be translated into an sbatch request like: `--time <time>`. Default: 0-05:59

* [`--partitions`] The partitions to request. This will be translated into an sbatch request like: `-p <partitions>`. default: general (N.B. This value should be a comma-separated list with no spaces, for example: `--partitions general,long`)

* [`--no-output`] If this flag is present, stdout will be redirected to /dev/null. This will be translated into an sbatch request like: `--output=/dev/null`. Default: If the flag is not present, then stdout will be directed to a log file with the job number. This corresponds to `--output=slurm-%J.out` in the sbatch call.

* [`--no-error`] If this flag is present, stderr will be redirected to /dev/null. This will be translated into an sbatch request like: `--error=/dev/null`. Default: If the flag is not present, then stderr will be directed to a log file with the job number. This corresponds to `--error=slurm-%J.err` in the sbatch call.

* [`--stdout-file`] If this is present and the --no-output flag is not given, then stdout will be directed to this file. This corresponds to `--output=stdout-file` in the sbatch call. Default: `slurm-%J.out`

* [`--stderr-file`] If this is present and the --no-error flag is not given, then stderr will be directed to this file. This corresponds to `--error=stderr-file` in the sbatch call. Default: `slurm-%J.err`

* [`--mail-type`] When to send an email notification of the job status. See official documentation for a description of the values. If a mail-user is not specified, this will revert to 'None'. Defaut: FAIL TIME_LIMIT

* [`--mail-user`] To whom an email will be sent, in accordance with mail-type. Default: None


Some example calls:

```
# This will submit the prepare-genome script to SLURM as a single job. That job
# will request 10 CPUs and 100G of RAM.

prepare-rpbp-genome config.yaml --num-cpus 10 --mem 100G --overwrite --logging-level INFO --use-slurm
```

```
# This will submit each sample as a separate job to SLURM. Each submitted job will request 10 cpus
# and 100G of RAM. 

# For example (refer to the example dataset), if c-elegans-test.yaml specifies 2 samples in the riboseq_samples
# value, then 2 jobs will be submitted to SLURM, one for each sample.

# If the --merge-replicates flag is given, then all of the profiles are first created. Then, the
# combined profiles are created in accordance with riboseq_biological_replicates from the config
# file. Finally, the last phase of the pipeline (predict-translated-orfs) is called.

# If any step fails, then later phases will not be started.

run-all-rpbp-instances c-elegans-test.yaml --num-cpus 10 --mem 100G --merge-replicates --use-slurm --logging-level INFO
```

[Back to top](#toc)