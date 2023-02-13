import sys

"""Define all default parameters and options.

** NOTE: Modifying the default parameter values set in this file
will have NO effect. To override a given parameter value, you must
either provide it via command argument when calling the Rp-Bp
pipeline, or else set the desired value in the configuration file.
Please refer to the documentation: https://rp-bp.readthedocs.io/en/latest/usage-instructions.html#


Default options for external programs (Flexbar, STAR) are
overridden via command line. Currently, call to Bowtie2 is not customisable.

Default parameters for Rp-Bp (default options for Rp-Bp scripts, including
shared MCMC options, but not processing options) are overridden by
providing the option key-value pair in the configuration file.

Processing options (parallel processing, logging options) are given via
command line.
"""


# default: processing
# overridden via command line arguments

default_num_cpus = 1
default_mem = "2G"

default_num_groups = 100  # currently cannot be overridden


# default: Rp-Bp (genome index creation, ORF identification)
# overridden via config file

default_start_codons = ["ATG"]
default_stop_codons = ["TAA", "TGA", "TAG"]


# default: Rp-Bp (metagene, ORF profiles)
# overridden via command line arguments

default_read_files_command = "zcat"
if sys.platform.startswith("darwin"):
    default_read_files_command = "gzcat"

# NOTES:
#   --limitBAMsortRAM 0 (with --genomeLoad NoSharedMemory)
#   limitBAMsortRAM is set to args.mem or default_mem, if not given, at run time

#   --outTmpDir is set to args.tmp at run time, else default STAR tmp is used (current directory)
#   if given, STAR args.tmp will be created

star_executable = "STAR"
# default options for `--runMode genomeGenerate`
star_options_genome = {}
# default options for mapping
star_options = {
    "readFilesCommand": default_read_files_command,
    "limitBAMsortRAM": 0,
    "alignIntronMin": 20,
    "alignIntronMax": 100000,
    "outFilterMismatchNmax": 1,
    "outFilterMismatchNoverLmax": 0.04,  # would require a minimum *mapped* length of 25 given max mismatch of 1
    "outFilterType": "BySJout",
    "outFilterIntronMotifs": "RemoveNoncanonicalUnannotated",
    "outSAMattributes": ["AS", "NH", "HI", "nM", "MD"],
    "outSAMtype": "BAM SortedByCoordinate",
    "sjdbOverhang": 33,  # roughly 90 percentile of a large dataset of varying Ribo-seq fragment lengths
    "seedSearchStartLmaxOverLread": 0.5,  # default seedSearchStartLmax normalised to read length
    "winAnchorMultimapNmax": 100,  # increase number of loci anchors are allowed to map to
}

# leave outFilterMultimapNmax to default 20, we filter the multimappers afterwards if desired

flexbar_options = {
    "max-uncalled": 1,
    "pre-trim-left": 0,
    "qtrim-format": "sanger",
    "qtrim": "TAIL",
    "qtrim-threshold": 10,
    "zip-output": "GZ",
}


# default: Rp-Bp (metagene, ORF profiles)
# overridden via config file

# Stan model instantiation
# Stan's multi-threaded processing is based on the Intel Threading Building Blocks (TBB) library,
# which must be linked to by the C++ compiler. True uses cpp_options={'STAN_THREADS': 'TRUE'}, or False.
model_inst_options = {"stan_threads": False}

# shared options
mcmc_shared = {
    "seed": 8675309,
    "chains": 2,
}

metagene_options = {
    "metagene_start_upstream": 300,
    "metagene_start_downstream": 300,
    "metagene_end_upstream": 300,
    "metagene_end_downstream": 300,
    "periodic_offset_start": -20,
    "periodic_offset_end": 0,
    "metagene_profile_length": 21,
    "seed": mcmc_shared["seed"],
    "chains": mcmc_shared["chains"],
    "metagene_iterations": 500,  # incl. warmup
    "min_metagene_profile_count": 1000,
    "min_metagene_image_count": 500,  # TODO: profiles with count < min_metagene_image_count will not be displayed (only for report)
    "min_metagene_bf_mean": 5,
    "max_metagene_bf_var": None,
    "min_metagene_bf_likelihood": 0.5,
}


# default: Rp-Bp (translation prediction)
# overridden via config file

translation_options = {
    "smoothing_fraction": 0.2,  # defaults are not written to the file names
    "smoothing_reweighting_iterations": 0,
    "orf_min_length_pre": 0,  # ORF with length < orf_min_length_pre are not processed, 0 ignore option
    "orf_max_length_pre": 0,  # ORF with length > orf_max_length_pre are not processed, 0 ignore option
    "orf_min_profile_count_pre": 5,  # ORF with profile sum < orf_min_profile_count_pre are not processed
    "seed": mcmc_shared["seed"],
    "chains": mcmc_shared["chains"],
    "translation_iterations": 500,  # incl. warmup
    "min_bf_mean": 5,
    "min_bf_likelihood": 0.5,
    "max_bf_var": None,
    "orf_min_length": 8,  # ORF > orf_min_length (nucleotides) are kept in the final set
    "orf_min_profile_count": None,  # ORF with sum across all frame > orf_min_profile_count are kept
    "chisq_alpha": 0.01,  # only relevant with the Chi2 pipeline
}


# default: Rp-Bp ORF labels

# Do NOT modify keys, they are used to label ORFs
orf_type_name_map = {
    "canonical": "CDS",
    "canonical_variant": "altCDS",
    "internal": "intORF",
    "five_prime": "uORF",
    "five_prime_overlap": "uoORF",
    "three_prime": "dORF",
    "three_prime_overlap": "doORF",
    "noncoding": "ncORF",
    "suspect": "Suspect",
    "overlap": "Overlap",
    "novel": "Novel",  # intergenic - not contained in the annotations
    "novel_canonical_variant": "Novel altCDS",
    "novel_noncoding": "Novel ncORF",  # is this a valid label?
    "novel_overlap": "Novel overlap",  # include any overlap e.g. 5'/3'
    "novel_suspect": "Novel suspect",
    "novel_internal": "Novel intORF",  # these 3 are in principle invalid labels
    "novel_five_prime": "Novel uORF",
    "novel_three_prime": "Novel dORF",
}

orf_type_labels = {
    "canonical": ["canonical"],
    "canonical_variant": ["canonical_variant", "internal", "novel_canonical_variant"],
    "five_prime": ["five_prime", "five_prime_overlap"],
    "three_prime": ["three_prime", "three_prime_overlap"],
    "noncoding": ["noncoding", "novel_noncoding"],
    "novel": ["novel"],
    "other": [
        "overlap",
        "novel_overlap",
        "suspect",
        "novel_suspect",
        "novel_internal",
        "novel_five_prime",
        "novel_three_prime",
    ],
}

# color blind hex colours
orf_type_colors = {
    "#0173b2": "canonical",
    "#56b4e9": "canonical_variant",
    "#029e73": "five_prime",
    "#de8f05": "three_prime",
    "#cc78bc": "noncoding",
    "#ece133": "novel",
    "#949494": "other",
}
