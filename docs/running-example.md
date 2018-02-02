
# Running Rp-Bp on the example dataset

<a id='toc'></a>

* [Download the dataset](#download)
* [Content of the file](#example-dataset-files)
* [Running Rp-Bp on the example](#running-example)
    * [Creating the reference index files](#creating-reference-indices)
    * [Running the main pipeline (with replicates)](#running-rpbp-pipeline)
    * [Troubleshooting](#common-problems)

---

<a name="download"></a>

## Download the dataset

A small example dataset using _C. elegans_ is available to [download](http://cloud.dieterichlab.org/index.php/s/7XHsCqZqU9AbQqB/download). Alternatively, the 
following commands can be used to download and extract the example .tar.gz file:

```
wget http://cloud.dieterichlab.org/index.php/s/7XHsCqZqU9AbQqB/download -O c-elegans-chrI-example.tar.gz
tar -xvf c-elegans-chrI-example.tar.gz
```

<a id="example-dataset-files"></a>

## Content of the file

The example dataset is distributed as a .tar.gz file and includes the following:

* `c-elegans-test.yaml`. The configuration file, used for creating the reference index files and for running the prediction pipeline. It includes all default options for creating the indices and for running the main pipeline, as well as detailed descriptions. 
The **exception** is the `min_metagene_profile_count` option, which has a value of 10 rather than its default of 1000. This is set artificially low because of the small number of reads in the sample dataset. Similarly, when plotting the results, `--min-visualization-count` has to be set to a lower value, see [Preprocessing analysis](analysis-script.md#preprocessing-report).
* `WBcel235.chrI.fa`. The reference sequence of Chromosome I for _C. elegans_.
* `WBcel235.79.chrI.gtf`. The Ensembl, version 79 annotations for Chromosome I for _C. elegans_.
* `X03680_1.fasta`. The sequences of the ribosomal subunits for _C. elegans_. The reference accession is X03680.1.
* `riboseq-adapters.fa`. An example adapter file to use with `flexbar`. It includes typical TruSeq and ArtSeq adapters, as well as a few adapters from the literature. It also includes a custom adapter used to create the sample dataset.
* `c-elegans.test-chrI.rep-1.fastq.gz`. A small test Ribo-seq dataset. It has been constructed to include some reads which uniquely map to the annotated transcripts, some reads which map to ribosomal sequences, some reads which do not uniquely map to the genome and some reads which are filtered due to quality issues.
* `c-elegans.test-chrI.rep-2.fastq.gz`. Another small test Ribo-seq dataset.
* `expected-orf-predictions`. The expected predictions and sequence files for each replicate (c-elegans-rep-1 and c-elegans-rep-2 files) and the merged replicates (c-elegans-test files). For an explanation of the output files and format, see [Predicting translated open reading frames](usage-instructions.md#predicting-translated-open-reading-frames).

Comparing your results with the expected output can be done *e.g.* using `bedtools` with the BED files containing the list of predicted ORFs. 

Due to differences among versions of the external programs used in the pipeline, it is however unlikely that all files will match exactly. If these differ significantly, it is possible that something is not working correctly in the pipeline. In such case, you can run the pipeline using the "DEBUG" logging level (see the [usage instructions](usage-instructions.md#logging-options)). This causes the scripts to output detailed runtime information which can be helpful for tracking down problems. If the problem is still not clear, please report the problem at the [github bug tracker](https://github.com/dieterich-lab/rp-bp/issues).

[Back to top](#toc)

<a id="running-example"></a>

## Running Rp-Bp on the example

**Before running the example** the paths in the `c-elegans-test.yaml` configuration file for each items below must be updated (*i.e.* `/path/to/your/c-elegans-example/` needs to be changed to point to the correct location/files).

* `genome_base_path`
* `gtf`
* `fasta`
* `ribosomal_fasta`
* `adapter_file`

* `ribosomal_index`
* `star_index`
* `riboseq_samples`
* `riboseq_data`

<a id='creating-reference-indices'></a>

### Creating the reference index files

The following command will create the necessary reference indices:
 
```
prepare-rpbp-genome c-elegans-test.yaml [--overwrite] [logging options] [processing options]
```

#### Command line options

* [`--overwrite`] Unless this flag is given, then steps for which the output files already exist will be skipped.
* [`logging options`] See [logging options](#logging-options).
* [`processing options`] See [parallel processing options](#parallel-processing-options).

**N.B The script reads all of the required paths from the configuration file, so it is important that all the paths point to the correct locations, as explained above.**

In total, creating the reference index files should take about 5 minutes.


[Back to top](#toc)

<a id='running-rpbp-pipeline'></a>

### Running the Rp-Bp pipeline

Example calls:

```
# Do not merge the replicates.
run-all-rpbp-instances c-elegans-test.yaml --overwrite --num-cpus 2 --logging-level INFO --keep-intermediate-files

# Merge the replicates, do not calculate Bayes factors nor make predictions for individual samples.
run-all-rpbp-instances c-elegans-test.yaml --overwrite --num-cpus 2 --logging-level INFO --merge-replicates --keep-intermediate-files

# Merge the replicates and also calculate Bayes factors and make predictions for individual samples.
run-all-rpbp-instances c-elegans-test.yaml --overwrite --num-cpus 2 --logging-level INFO --merge-replicates --run-replicates --keep-intermediate-files
```

The call to `run-all-rpbp-instances` requires the configuration file. All other arguments are optional. In all the example calls above, the command will run the Rp-Bp complete translation prediction pipeline using 2 CPUS. For more details regarding the options and the expected output files, please consult the [usage instructions](usage-instructions.md#running-pipelines).

N.B. The `--overwrite` flag is given above to ensure that all of the files are (re-)created. In typical use cases, if some of the files already exist (*e.g.* the quality-filtered reads), then this flag can be omitted.

**Using replicates**

The Rp-Bp pipeline handles replicates by adding the (smoothed) ORF profiles. The Bayes factors and predictions are then calculated based on the combined profiles. The `--merge-replicates` flag indicates that the replicates should be merged. By default, if the `--merge-replicates` flag is given, then predictions will not be made for the individual datasets, unless the `--run-replicates` flag is also given, in which case predictions will be made for both the merged replicates as well as the individual samples.

The replicates are specified by `riboseq_biological_replicates` in the configuration file. This value should be a dictionary, where the key of the dictionary is a string description of the condition and the value is a list that gives all of the sample replicates which belong to that condition. The names of the sample replicates must match the sample names specified in `riboseq_samples`. 


[Back to top](#toc)

<a id='common-problems'></a>

### Troubleshooting

Some common problems result due to versions of external programs. The can be controlled using command line options to `run-all-rpbp-instances`.


* `--star-executable`. In principle, `STARlong` (as opposed to `STAR`) could be used for alignment. Given the nature of riboseq reads (that is, short due to the experimental protocols of degrading everything not protected by a ribosome), this is unlikely to be a good choice, though. Default: `STAR`


* `--star-read-files-command`. The input for `STAR` will always be a gzipped fastq file. `STAR` needs the system command which means "read a gzipped text file". As discovered in [Issue #35](https://github.com/dieterich-lab/rp-bp/issues/35), the name of this command is different on OSX and Ubuntu. The program now attempts to guess the name of this command based on the operating system, but it can be explicitly specified as a command line option. Default: `gzcat` if `sys.platform.startswith("darwin")`; `zcat` otherwise. Please see [python.sys documentation](https://docs.python.org/3/library/sys.html) for more details about attempting to guess the operating system.

[Back to top](#toc)
