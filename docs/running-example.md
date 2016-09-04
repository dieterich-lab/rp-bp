
# Running the small example dataset

A small example dataset using _C. elegans_ is available for [download](http://cloud.dieterichlab.org/index.php/s/3cyluM3ZCsvf0PT/download). Please see [below](#example-dataset-files) for the exact contents of the download, as well as instructions for downloading it from the command line.

Additionally, the expected outputs of the pipeline are included. Due to differences among versions of the external programs used in the pipeline (samtools, etc.), it is unlikely that all intermediate files will match exactly. However, we do include a script to compare the ORFs predicted as translated using the pipeline to those which are expected. If these differ significantly, it suggests something is not working correctly in the pipeline.

If the results differ significantly, please run the pipeline using the "DEBUG" logging level (see the [usage instructions](usage-instructions.ipynb#logging-options)). This causes the scripts to output detailed runtime information which can be helpful for tracking down problems. If the problem is still not clear, please report the problem at the [github bug tracker](https://github.com/dieterich-lab/rp-bp/issues).

In total, creating the reference index files should take about 5 minutes and running the main pipeline should take an additional 15 to 20 minutes on a commodity laptop.

<a id='toc'</a>
* [Example dataset files](#example-dataset-files)
* [Creating the reference index files](#creating-reference-indices)
* [Running the Rp-Bp pipeline (also with replicates)](#running-rpbp-pipeline)
* [Validating the Rp-Bp prediction results](#validating-results)

<a id="example-dataset-files"></a>

## Example dataset files

The example dataset is distributed as a .tar.gz file and includes the following:
* `WBcel235.79.chrI.yaml`. The configuration file for creating the reference index files. It includes all possible options for creating the indices as well as detailed descriptions.

* `WBcel235.dna.toplevel.chrI.fa`. The reference sequence of Chromosome I for _C. elegans_.

* `WBcel235.79.chrI.gtf`. The Ensembl, version 79 annotations for Chromosome I for _C. elegans_.

* `X03680_1.fasta`. The sequences of the ribosomal subunits for _C. elegans_. The reference accession is X03680.1.

* `c-elegans-test.yaml`. The configuration file for running the prediction pipeline. This example configuration file includes all possible options for the pipeline with detailed explanations of the options. The **exception** is the `min_metagene_profile_count` option, which has a value of 10 rather than its default of 1000. This is set artificially low because of the small number of reads in the sample dataset.

* `riboseq-adapters.fa`. An example adapter file for use with `flexbar`. It includes typical TruSeq and ArtSeq adapters, as well as a few adapters from the literature. It also includes a custom adapter used to create the sample dataset.

* `c-elegans.test-chrI.rep-1.fastq.gz`. A small test sequencing dataset. It has been constructed to include some reads which uniquely map to the annotated transcripts, some reads which map to ribosomal sequences, some reads which do not uniquely map to the genome and some reads which are filtered due to quality issues.

* `c-elegans.test-chrI.rep-2.fastq.gz`. Another small test sequencing dataset.

* `expected-orf-predictions`. The expected predictions and sequence files for each replicate (c-elegans-rep-1 and c-elegans-rep-2 files) and the merged replicates (c-elegans-test files). Please see the [usage instructions](usage-instructions.ipynb#logging-options) for the meaning of each of the files.

**Downloading from the command line**

The following `wget` command can be used to download the example .tar.gz file:


```python
wget http://cloud.dieterichlab.org/index.php/s/3cyluM3ZCsvf0PT/download -O c-elegans-chrI-example.tar.gz
```

[Back to top](#toc)

<a id='creating-reference-indices'></a>

## Creating the reference index files

**Before running the example** the paths in the `WBcel235.79.yaml` configuration file must be updated to point to the correct locations. The following configuration values should be updated to point to the appropriate files in the example. (Mostly, `/home/bmalone/python-projects/rp-bp/data/` should be replaced to the location of the examples.)

* `gtf`
* `fasta`
* `ribosomal_fasta`
* `genome_base_path`
* `ribosomal_index`
* `star_index`

The following command will create the necessary reference files using 2 CPUS and 4GB of RAM for STAR. Please see the [usage instructions](usage-instructions.ipynb#creating-reference-genome-indices) for the expected output files.

The `--use-slurm` and related options can also be used if SLURM is available. Please see the [usage instructions](usage-instructions.ipynb#parallel-processing-options) for more information.

N.B. The `--overwrite` flag is given below to ensure all of the files are (re-)created. In typical use cases, if some of the files already exist (e.g., the STAR index), then this flag can be omitted.

This command should only take about 5 minutes on recent commodity hardware (such as a laptop).


```python
prepare-genome WBcel235.79.chrI.yaml --num-cpus 2 --mem 4G --overwrite --logging-level INFO
```

[Back to top](#toc)

<a id='running-rpbp-pipeline'></a>

## Running the Rp-Bp pipeline

**Before running the example** the paths in the `c-elegans-test.yaml` configuration file must be updated to point to the correct locations. The following configuration values should be updated to point to the appropriate files in the example. (Mostly, `/home/bmalone/python-projects/rp-bp/data/` should be replaced to the location of the examples.)

Reference files and locations should be exactly the same as used in the  `WBcel235.79.yaml` file.

* `gtf`
* `fasta`
* `genome_base_path`
* `ribosomal_index`
* `star_index`

The sample and output file paths must also be updated.

* `riboseq_samples`
* `riboseq_data`
* `adapter_file`

The following command will run the Rp-Bp (and Rp-chi) translation prediction pipelines using 2 CPUS. Please see the [usage instructions](usage-instructions.ipynb#running-pipelines) for the expected output files.

The `--use-slurm` and related options can also be used if SLURM is available. Please see the [usage instructions](usage-instructions.ipynb#parallel-processing-options) for more information.

N.B. The `--overwrite` flag is given below to ensure all of the files are (re-)created. In typical use cases, if some of the files already exist (e.g., the quality-filtered reads), then this flag can be omitted.

N.B. While performing the MCMC sampling, many messages indicating the "Elapsed Time" will be printed. This is a [known issue](https://github.com/stan-dev/pystan/issues/98) with pystan. Additionally, many "Informational Message: The current Metropolis proposal is about to be rejected because of the following issue" may also appear. These are also expected and (typically) do not indicate an actual problem.

**Using replicates**

The Rp-Bp pipeline handles replicates by adding the (smoothed) ORF profiles. The Bayes factors and predictions are then calculated based on the combined profiles. The `--merge-replicates` flag indicates that the replicates should be merged. By default, if the `--merge-replicates` flag is given, then predictions will not be made for the individual datasets. The `--run-replicates` flag can be given to override this and make predictions for both the merged replicates as well as the individual datasets.

The replicates are specified by `riboseq_biological_replicates` in the configuration file. This value should be a dictionary, where the key of the dictionary is a string description of the condition and the value is a list that gives all of the sample replicates which belong to that condition. The names of the sample replicates must match the dataset names specified in `riboseq_samples`.


```python
# running on each dataset separately
process-all-samples c-elegans-test.yaml --overwrite --num-cpus 2 --logging-level INFO

# merging the replicates, do not calculate Bayes factors and make predictions for individual datasets
process-all-samples c-elegans-test.yaml --overwrite --num-cpus 2 --logging-level INFO --merge-replicates

# merging the replicates and calculating Bayes factors and making predictions for individual datasets
process-all-samples c-elegans-test.yaml --overwrite --num-cpus 2 --logging-level INFO --merge-replicates --run-replicates
```

[Back to top](#toc)

<a id='validating-results'></a>

## Validating the Rp-Bp prediction results

As described above, we do not necessarily expect every intermediate file to exactly match across different versions of the external programs. However, the final predictions should be rather similar across versions.

With this in mind, we include a script with the Rp-Bp package which calculates the difference, in base pairs, among the ORFs predicted as translated between the "expected" outputs for the example dataset and those calculated on a particular run. While no hard threshold distinguishes between "successful" runs or not, the value should not be very large.

The following command will print the number of unique bases to the predictions and expected predictions, as well as the overlap between them. The filenames may need to be changed to account for the local paths.


```python
# to check the predictions from Rp-Bp
calculate-bed-overlap c-elegans-test.expected.rpbp.predicted-orfs.bed.gz orf-predictions/c-elegans-chrI.test-unique.length-29.offset-12.predicted-orfs.bed.gz

# to check the predictions from Rp-chi
calculate-bed-overlap c-elegans-test.expected.rpchi.predicted-orfs.bed.gz orf-predictions/c-elegans-chrI.test-unique.length-29.offset-12.chisq.predicted-orfs.bed.gz
```

[Back to top](#toc)
