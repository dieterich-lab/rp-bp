
import sys

"""Define all default parameters and options.

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
default_mem = '2G'

default_num_groups = 100 # currently cannot be overridden


# default: Rp-Bp (genome index creation, ORF identification)
# overridden via config file

default_start_codons = ['ATG']
default_stop_codons = ['TAA', 'TGA', 'TAG']


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

star_executable = 'STAR'
star_options = {
    'readFilesCommand': default_read_files_command,
    'limitBAMsortRAM': 0,
    'alignIntronMin': 20,
    'alignIntronMax': 100000,
    'outFilterMismatchNmax': 1,
    'outFilterMismatchNoverLmax': 0.04,
    'outFilterType': 'BySJout',
    'outFilterIntronMotifs': 'RemoveNoncanonicalUnannotated',
    'outSAMattributes': ['AS', 'NH', 'HI', 'nM', 'MD'],
    'outSAMtype': 'BAM SortedByCoordinate',
    'sjdbOverhang': 50
}

flexbar_options = {
    'max-uncalled': 1,
    'pre-trim-left': 0,
    'qtrim-format': 'sanger',
    'qtrim': 'TAIL',
    'qtrim-threshold': 10,
    'zip-output': 'GZ'
}


# default: Rp-Bp (metagene, ORF profiles)
# overridden via config file

metagene_options = {
    'start-upstream': 300,
    'start-downstream': 300,
    'end-upstream': 300,
    'end-downstream': 300,
    'periodic-offset-start': -20,
    'periodic-offset-end': 0,
    'metagene-profile-length': 21,
    'seed': 8675309,
    'chains': 2,
    'iterations': 500
}


# default: Rp-Bp (translation prediction)
# overridden via config file

translation_options = {
    'fraction': 0.2,  # defaults are not written in the filenames
    'reweighting-iterations': 0,
    'min-length-pre': 0, # ORF with length < min-length are not processed, 0 ignore option
    'max-length-pre': 0, # ORF with length > max-length are not processed, 0 ignore option
    'min-profile-pre': 5, # ORF with profile sum < min-profile are not processed
    'seed': 8675309,
    'chains': 2,
    'iterations': 500,
    'orf-types': [], # predict only these, if empty predict all types
    'min-bf-mean': 5,
    'min-bf-likelihood': 0.5,
    'max-bf-var': None,
    'min-length': 20,  # ORF > min-length (nucleotides) are kept in the final set
    'min-profile': None,    # ORF with sum across all frame > minimum-profile-sum are kept
    'chisq-sig': 0.01 # only relevant with the Chi2 pipeline
}
