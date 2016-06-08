import glob
import os

### parameterized names

def get_cds_only_string(is_cds_only):
    cds_only = ""
    if is_cds_only:
        cds_only = ".cds-only"
    return cds_only

def get_chisq_string(is_chisq):
    chisq = ""
    if is_chisq:
        chisq = ".chisq"
    return chisq

def get_fastqc_name(filename):
    """ Given the sequence or alignment filename, this function extracts the name
        used by fastqc. In particular, it removes the ".fasta[.gz]" or ".fastq[.gz]"
        or ".sam" or ".bam" ending from the filename.

        Args:
            filename: the name of the sequence or alignment file (NOT including the path)

        Returns:
            string : the fastqc name

        Imports:
            misc.utils
    """

    import misc.utils as utils

    # first, get the filename
    filename = utils.get_basename(filename)
    filename = filename.replace(".fasta", "")
    filename = filename.replace(".fastq", "")
    filename = filename.replace(".sam", "")
    filename = filename.replace(".bam", "")
    filename = filename.replace(".gz", "")

    return filename


def get_length_string(length=None):
    l = ""
    if length is not None:
        if isinstance(length, (list, tuple)):
            l = "-".join(length)
            l = ".length-{}".format(l)
        else:
            l = ".length-{}".format(length)
    return l

def get_merged_string(is_merged):
    m = ""
    if is_merged:
        m = ".merged"
    return m

def get_micro_string(is_micro):
    m = ""
    if is_micro:
        m = ".micro-only"
    return m

def get_note_string(note=None):
    note_str = ""
    if (note is not None) and  (len(note) > 0):
        note_str = ".{}".format(note)
    return note_str

def get_offset_string(offset=None):
    o = ""
    if offset is not None:
        if isinstance(offset, (list, tuple)):
            o = "-".join(offset)
            o = ".offset-{}".format(o)
        else:
            o = ".offset-{}".format(offset)
    return o


def get_unique_string(is_unique):
    unique = ""
    if is_unique:
        unique = "-unique"
    return unique

def get_transcriptome_string(is_transcriptome):
    transcriptome = ""
    if is_transcriptome:
        transcriptome = ".transcriptome"
    return transcriptome


### m

# used
def get_metagene_profiles(riboseq_base, name, length=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, note=None):

    s = get_riboseq_base(riboseq_base, name, 'metagene-profiles', length=length, 
            is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, note=note)
    s = s + ".metagene-profile.csv.gz"
    return s

# used
def get_metagene_profiles_bayes_factors(riboseq_base, name, length=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, note=None):

    s = get_riboseq_base(riboseq_base, name, 'metagene-profiles', length=length, 
            is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, note=note)
    s = s + ".metagene-periodicity-bayes-factors.csv.gz"
    return s


# used
def get_models(models_base, model_type):
    path_ex = os.path.join(models_base, model_type, '*pkl')
    models = glob.glob(path_ex)
    return models

### o

# used
def get_orfs(base_path, name, note=None):
    note_str = get_note_string(note)
    fn = '{}.genomic-orfs{}.bed.gz'.format(name, note_str)
    return os.path.join(base_path, 'transcript-index', fn)

### p

# used
def get_peptide_coverage_line_graph(riboseq_base, name, length=None, offset=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, note=None, image_type='pdf'):
    
    subfolder = os.path.join('peptide-matches', 'plots')
    s = get_riboseq_base(riboseq_base, name, subfolder, length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, 
        note=note)
    s = s + ".orf-peptide-coverage.{}".format(image_type)
    return s


# used
def get_periodic_offsets(riboseq_base, name, is_unique=False, is_cds_only=False, is_transcriptome=False, note=None):
    s = get_riboseq_base(riboseq_base, name, 'metagene-profiles', is_unique=is_unique, 
        is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, note=note)

    s = s + ".periodic-offsets.csv.gz"
    return s

def get_preprocessing_report(base_path):
    fn = "preprocessing-report.tex"
    return os.path.join(base_path, 'preprocessing-report', fn)

### r
def get_raw_data_path(base_path):
    return os.path.join(base_path, 'raw-data')

def get_raw_data_fastqc_path(base_path):
    rdp = get_raw_data_path(base_path)
    return os.path.join(rdp, 'fastqc')


def get_raw_data_fastqc_data(base_path, filename):

    name = get_fastqc_name(filename)
    fastqc_folder = '{}_fastqc'.format(name)
    rdp = get_raw_data_fastqc_path(base_path)

    p = os.path.join(rdp, fastqc_folder, 'fastqc_data.txt')
    #print("raw data path: {}".format(rdp))
    #print("fastqc folder: {}".format(fastqc_folder))
    #print("hey!!! i'm here!!! p: {}".format(p))
    return p

### riboseq

# b

# used
def get_riboseq_base(riboseq_base, name, sub_folder, length=None, offset=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, is_chisq=False, note=None):
    
    cds_only = get_cds_only_string(is_cds_only)
    unique = get_unique_string(is_unique)
    l = get_length_string(length)
    o = get_offset_string(offset)
    transcriptome = get_transcriptome_string(is_transcriptome)
    chisq = get_chisq_string(is_chisq)
    n = get_note_string(note)
    return os.path.join(riboseq_base, sub_folder, 
            '{}{}{}{}{}{}{}{}'.format(name, n, transcriptome, unique, cds_only, l, o, chisq))

# used
def get_riboseq_bam_base(riboseq_base, name, length=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, note=None):

    bam_base = get_riboseq_base(riboseq_base, name, 'without-rrna-mapping', length=length, is_unique=is_unique, 
        is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, note=note)
    return bam_base

# used
def get_riboseq_bam(riboseq_base, name, length=None, is_unique=False, is_cds_only=False, 
        is_transcriptome=False, note=None):

    s = get_riboseq_bam_base(riboseq_base, name, length=length, is_unique=is_unique, 
            is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, note=note)
    s = s + ".bam"
    return s

def get_riboseq_bam_fastqc_path(riboseq_data):
    return os.path.join(riboseq_data, 'without-rrna-mapping', 'fastqc')

def get_riboseq_bam_fastqc_data(riboseq_data, name, length=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, note=None):

    cds_only = get_cds_only_string(is_cds_only)
    unique = get_unique_string(is_unique)
    l = get_length_string(length)
    transcriptome = get_transcriptome_string(is_transcriptome)
    n = get_note_string(note)
    name = '{}{}{}{}{}{}'.format(name, n, transcriptome, unique, cds_only, l)

    fastqc_folder = '{}_fastqc'.format(name)
    return os.path.join(riboseq_data, 'without-rrna-mapping', 'fastqc', fastqc_folder, 'fastqc_data.txt')
    
# used
def get_riboseq_bayes_factors(riboseq_base, name, length=None, offset=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, note=None):
    
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, 
            is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, note=note)
    
    s = s + ".bayes-factors.bed.gz"
    return s

# f
def get_riboseq_fastq(riboseq_data, name):
    return os.path.join(riboseq_data, 'raw-data', '{}.fastq.gz'.format(name))

# m

def get_riboseq_metagene_profile_image(riboseq_base, name, image_type='eps', length=None, is_unique=False, is_cds_only=False, is_transcriptome=False, note=None):
    s = get_riboseq_base(riboseq_base, name, 'metagene-profiles', length=length, is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, note=note)
    s = s + "." + image_type
    return s

def get_metagene_profile_bayes_factor_image(riboseq_base, name, image_type='eps', length=None, is_unique=False, is_cds_only=False, is_transcriptome=False, note=None):
    s = get_riboseq_base(riboseq_base, name, 'metagene-profiles', length=length, is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, note=note)
    s = s + ".bayes-factors." + image_type
    return s



# p

# used
def get_riboseq_peptide_matches(riboseq_base, name, length=None, offset=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, is_chisq=False, note=None):
    
    s = get_riboseq_base(riboseq_base, name, 'peptide-matches', length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, 
        is_chisq=is_chisq, note=note)
    s = s + ".peptide-matches.csv.gz"
    return s


# used
def get_riboseq_predicted_orfs(riboseq_base, name, length=None, offset=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, note=None, is_chisq=False):
    
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, 
        is_chisq=is_chisq, note=note)
    s = s + ".predicted-orfs.bed.gz"
    return s

# used
def get_riboseq_predicted_orfs_dna(riboseq_base, name, length=None, offset=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, is_chisq=False, note=None):
    
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, 
        is_chisq=is_chisq, note=note)
    s = s + ".predicted-orfs.dna.fa"
    return s

# used
def get_riboseq_predicted_orfs_protein(riboseq_base, name, length=None, offset=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, is_chisq=False, note=None):
    
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, 
        is_chisq=is_chisq, note=note)
    s = s + ".predicted-orfs.protein.fa"
    return s

def get_riboseq_predicted_orf_mackowiak_overlap_image(riboseq_base, name, image_type='eps', 
        length=None, offset=None, is_unique=False, is_cds_only=False, is_transcriptome=False, 
        is_chisq=False, note=None):
    
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, 
        is_chisq=is_chisq, note=note)
    s = s + ".predicted-orf-mackowiak-overlap.{}".format(image_type)
    return s

def get_riboseq_predicted_orf_peptide_coverage(riboseq_base, name, length=None, offset=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, note=None):
    
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, note=note)
    s = s + ".predicted-orf-peptide-coverage.csv.gz"
    return s


def get_riboseq_predicted_orf_peptide_coverage_image(riboseq_base, name, image_type='eps', length=None, 
        offset=None, is_unique=False, is_cds_only=False, is_transcriptome=False, note=None):
    
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, note=note)
    s = s + ".predicted-orf-peptide-coverage.{}".format(image_type)
    return s

def get_riboseq_predicted_orf_qti_seq_overlap_image(riboseq_base, name, image_type='eps', length=None, 
        offset=None, is_unique=False, is_cds_only=False, is_transcriptome=False, note=None):
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, note=note)
    s = s + ".predicted-orf-qti-seq-overlap.{}".format(image_type)
    return s


def get_riboseq_predicted_orf_spikes(riboseq_base, name, length=None, 
        offset=None, is_unique=False, is_cds_only=False, is_transcriptome=False, note=None):
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, note=note)
    s = s + ".spike-codons.csv.gz"
    return s
    
def get_riboseq_predicted_orf_spikes_bed(riboseq_base, name, length=None, 
        offset=None, is_unique=False, is_cds_only=False, is_transcriptome=False, note=None):
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, note=note)
    s = s + ".spikes.bed.gz"
    return s

def get_riboseq_predicted_orf_spikes_image(riboseq_base, name, image_type='eps', length=None, 
        offset=None, is_unique=False, is_cds_only=False, is_transcriptome=False, note=None):
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, note=note)
    s = s + ".spike-codons.{}".format(image_type)
    return s



def get_riboseq_predicted_orf_type_overlap_image(riboseq_base, name, image_type='eps', length=None, 
        offset=None, is_unique=False, is_cds_only=False, is_transcriptome=False, note=None):
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, note=note)
    s = s + ".predicted-orf-type-overlap.{}".format(image_type)
    return s


# used
def get_riboseq_profiles(riboseq_base, name, length=None, offset=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, note=None):

    s = get_riboseq_base(riboseq_base, name, 'orf-profiles', length=length, offset=offset, 
            is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, note=note)

    s = s + ".profiles.mtx"
    return s

# r

def get_riboseq_read_filtering_counts(riboseq_base):
    s = os.path.join(riboseq_base, 'read-filtering-counts.csv.gz')
    return s

def get_riboseq_read_filtering_counts_image(riboseq_base, note="", image_type="eps"):
    note_str = get_note_string(note)
    fn = "read-filtering-counts.{}{}".format(note_str, image_type)
    s = os.path.join(riboseq_base, fn)
    return s

    

### s

# used
def get_star_index(base_path, name, is_merged=False):
    m = get_merged_string(is_merged)
    fn = '{}{}'.format(name, m)
    return os.path.join(base_path, 'STAR', fn)


### t

# used
def get_transcript_fasta(base_path, name):
    fn = '{}.transcripts.fa'.format(name)
    return os.path.join(base_path, 'transcript-index', fn)


### w

# used
def get_without_adapters_base(base_path, name, note=None):
    n = get_note_string(note)
    base = "{}{}".format(name, n)
    return os.path.join(base_path, 'without-adapters', base)

# used
def get_without_adapters_fastq(base_path, name, note=None):
    n = get_note_string(note)
    return os.path.join(base_path, 'without-adapters', '{}{}.fastq.gz'.format(name, n))

def get_without_adapters_fastqc(base_path):
    return os.path.join(base_path, 'without-adapters', 'fastqc')

def get_without_adapters_fastqc_data(base_path, name, note=None):
    n = get_note_string(note)
    fastqc_folder = '{}{}_fastqc'.format(name, n)
    return os.path.join(base_path, 'without-adapters', 'fastqc', fastqc_folder, 'fastqc_data.txt')


# used
def get_with_rrna_fastq(base_path, name, note=None):
    n = get_note_string(note)
    return os.path.join(base_path, 'with-rrna', '{}{}.fastq.gz'.format(name, n))

# used
def get_with_rrna_fastqc(base_path):
    return os.path.join(base_path, 'with-rrna', 'fastqc')

def get_with_rrna_fastqc_data(base_path, name, note=None):
    n = get_note_string(note)
    fastqc_folder = '{}{}_fastqc'.format(name, n)
    return os.path.join(base_path, 'with-rrna', 'fastqc', fastqc_folder, 'fastqc_data.txt')

# used
def get_without_rrna_fastq(base_path, name, note=None):
    n = get_note_string(note)
    return os.path.join(base_path, 'without-rrna', '{}{}.fastq.gz'.format(name, n))

def get_without_rrna_fastqc(base_path):
    return os.path.join(base_path, 'without-rrna', 'fastqc')

def get_without_rrna_fastqc_data(base_path, name, note=None):
    n = get_note_string(note)
    fastqc_folder = '{}{}_fastqc'.format(name, n)
    return os.path.join(base_path, 'without-rrna', 'fastqc', fastqc_folder, 'fastqc_data.txt')

