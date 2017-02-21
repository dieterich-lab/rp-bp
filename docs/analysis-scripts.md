# Downstream analysis of the Rp-Bp results

## Creating read length-specific profiles

As described in the [usage instructions](usage-instructions.md#output-files-1), 
Rp-Bp writes the unsmoothed ORF profiles to a matrix market file. This profile 
merges reads of all lengths.

The `create-read-length-orf-profiles` script can be used to create profile files
which also include counts of individual read lengths.

```
create-read-length-orf-profiles <config> <sample or condition name> <out> [--is-condition]
```

### Command line options

* `config`. A yaml config file

* `sample or condition name`. The name of either one of the `riboseq_samples` or
  `riboseq_biological_replicates` from the config file

* `out`. The output (txt.gz) file, containing the read-length specific profiles.
  The format is a sparse coordinate format inspired by the [matrix market 
  format](http://math.nist.gov/MatrixMarket/formats.html). See below for 
  details about the output format.

* [`--is-condition`]. If the `sample or condition name` is a condition, that is,
  if it is a key from `riboseq_biological_replicates`, then this flag must be
  given.

Additionally, the command can be given [logging](usage-instructions.md#logging-options)
and [slurm](usage-instructions.md#parallel-processing-options) options.

### Output format

Each line in the output file is a tuple containing the following values.

* `read_length`. The (trimmed) read lengths for this position.

* `orf_num`. An identifier which maps to `orf_num` in the [(static) list of 
  ORFs](usage-instructions.md#output-files) for the reference, 
  `<genome_base_path>/transcript-index/<genome_name>.genomic-orfs.<orf_note>.bed.gz`.

* `orf_position`. The base-0 position with respect to the spliced transcript 
  (so `position % 3 == 0` implies the position is in-frame)

* `read_count`. The sum of counts across all replicates for the condition (if
  `--is-condition` is given) or the single sample (otherwise) after adjusting
    according to P-sites and removing multimappers.