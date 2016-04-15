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

def get_length_string(length=None):
    l = ""
    if length is not None:
        if isinstance(length, (list, tuple)):
            l = "-".join(length)
            l = ".length-{}".format(l)
        else:
            l = ".length-{}".format(length)
    return l

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





### c
def get_cds_bed(base_path, name):
    fn = '{}.transcript-cds-coordinates.bed'.format(name)
    return os.path.join(base_path, 'transcript-index', fn)

def get_cds_only_fasta(base_path, name):
    fn = '{}.transcripts.cds-only.fa'.format(name)
    return os.path.join(base_path, 'transcript-index', fn)

def get_cds_only_bowtie_index(base_path, name):
    fn = '{}.transcripts.cds-only'.format(name)
    return os.path.join(base_path, 'transcript-index', fn)

### g
def get_gtf_csv(base_path, name):
    fn = '{}.csv.gz'.format(name)
    return os.path.join(base_path, fn)

### o

def get_orfs(base_path, name, note=None):
    note_str = get_note_string(note)
    fn = '{}.genomic-orfs{}.bed.gz'.format(name, note_str)
    return os.path.join(base_path, 'transcript-index', fn)

### p
def get_preprocessing_report(base_path):
    fn = "preprocessing-report.tex"
    return os.path.join(base_path, 'preprocessing-report', fn)

### r
def get_raw_data_path(base_path):
    return os.path.join(base_path, 'raw-data')

def get_raw_data_fastqc_path(base_path):
    rdp = get_raw_data_path(base_path)
    return os.path.join(rdp, 'fastqc')


def get_raw_data_fastqc_data(base_path, name):
    fastqc_folder = '{}_fastqc'.format(name)
    rdp = get_raw_data_path(base_path)
    return os.path.join(rdp, 'fastqc', fastqc_folder, 'fastqc_data.txt')

### riboseq

# b
def get_riboseq_base(riboseq_base, name, sub_folder, length=None, offset=None, is_unique=False, is_cds_only=False, is_transcriptome=False, is_chisq=False, note=None):
    cds_only = get_cds_only_string(is_cds_only)
    unique = get_unique_string(is_unique)
    l = get_length_string(length)
    o = get_offset_string(offset)
    transcriptome = get_transcriptome_string(is_transcriptome)
    chisq = get_chisq_string(is_chisq)
    n = get_note_string(note)
    return os.path.join(riboseq_base, 'without-rrna-mapping', sub_folder, '{}{}{}{}{}{}{}{}'.format(name, n, transcriptome, unique, cds_only, l, o, chisq))

def get_riboseq_bam_base(riboseq_base, name, length=None, is_unique=False, is_cds_only=False, is_transcriptome=False):
    bam_base = get_riboseq_base(riboseq_base, name, 'bam', length=length, is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome)
    return bam_base

def get_riboseq_bam(riboseq_base, name, length=None, is_unique=False, is_cds_only=False, is_transcriptome=False, note=None):
    s = get_riboseq_bam_base(riboseq_base, name, length=length, is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome)
    s = s + ".bam"
    return s

def get_riboseq_bam_path(riboseq_data):
    return os.path.join(riboseq_data, 'without-rrna-mapping', 'bam')

def get_riboseq_bam_fastqc_path(riboseq_data):
    return os.path.join(riboseq_data, 'without-rrna-mapping', 'bam', 'fastqc')

def get_riboseq_bam_fastqc_data(riboseq_data, name, length=None, is_unique=False, is_cds_only=False, is_transcriptome=False, note=None):
    cds_only = get_cds_only_string(is_cds_only)
    unique = get_unique_string(is_unique)
    l = get_length_string(length)
    transcriptome = get_transcriptome_string(is_transcriptome)
    n = get_note_string(note)
    name = '{}{}{}{}{}{}'.format(name, n, transcriptome, unique, cds_only, l)

    fastqc_folder = '{}_fastqc'.format(name)
    return os.path.join(riboseq_data, 'without-rrna-mapping', 'bam', 'fastqc', fastqc_folder, 'fastqc_data.txt')
    
def get_riboseq_bayes_factors(riboseq_base, name, length=None, offset=None, is_unique=False, is_cds_only=False, is_transcriptome=False, note=None):
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, note=note)
    s = s + ".bayes-factors.bed.gz"
    return s


def get_riboseq_bitseq(riboseq_base, name, length=None, is_unique=False, is_cds_only=False, is_transcriptome=False):
    unique = get_unique_string(is_unique)
    cds_only = get_cds_only_string(is_cds_only)
    l = get_length_string(length)
    transcriptome = get_transcriptome_string(is_transcriptome)
    return os.path.join(riboseq_base, 'without-rrna-mapping', 'transcript-abundance', 'bitseq', '{}{}{}{}{}.bitseq'.format(name, transcriptome, unique, cds_only, l))

def get_riboseq_bitseq_malphas(riboseq_base, name, length=None, is_unique=False, is_cds_only=False, is_transcriptome=False):
    s = get_riboseq_bitseq(riboseq_base, name, length=length, is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome)
    s = s + ".m_alphas"
    return s


# f
def get_riboseq_fastq(riboseq_data, name):
    return os.path.join(riboseq_data, 'raw-data', '{}.fastq.gz'.format(name))

# g

def get_riboseq_gtf_signals(riboseq_base, name, length=None, is_unique=False, is_cds_only=False, is_transcriptome=False):
    s = get_riboseq_bam_base(riboseq_base, name, length=length, is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome)
    s = s + ".gtf.signals.csv.gz"
    return s

# m

def get_metagene_profiles(riboseq_base, name, length=None, is_unique=False, is_cds_only=False, is_transcriptome=False):
    s = get_riboseq_base(riboseq_base, name, 'metagene-profiles', length=length, is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome)
    s = s + ".metagene-profile.csv.gz"
    return s

def get_riboseq_metagene_profile_image(riboseq_base, name, image_type='eps', length=None, is_unique=False, is_cds_only=False, is_transcriptome=False):
    s = get_riboseq_base(riboseq_base, name, 'metagene-profiles', length=length, is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome)
    s = s + "." + image_type
    return s

def get_metagene_profile_bayes_factor_image(riboseq_base, name, image_type='eps', length=None, is_unique=False, is_cds_only=False, is_transcriptome=False):
    s = get_riboseq_base(riboseq_base, name, 'metagene-profiles', length=length, is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome)
    s = s + ".bayes-factors." + image_type
    return s


def get_metagene_profiles_bayes_factors(riboseq_base, name, length=None, is_unique=False, is_cds_only=False, is_transcriptome=False):
    s = get_riboseq_base(riboseq_base, name, 'metagene-profiles', length=length, is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome)
    s = s + ".metagene-periodicity-bayes-factors.csv.gz"
    return s


# p
def get_periodic_offsets(riboseq_base, name, is_unique=False, is_cds_only=False, is_transcriptome=False):
    s = get_riboseq_base(riboseq_base, name, 'metagene-profiles', is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome)
    s = s + ".periodic-offsets.csv.gz"
    return s

def get_riboseq_predicted_orfs(riboseq_base, name, length=None, offset=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, note=None, is_chisq=False):
    
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, is_chisq=is_chisq, note=note)
    s = s + ".predicted-orfs.bed.gz"
    return s

def get_riboseq_predicted_orfs_dna(riboseq_base, name, length=None, offset=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, is_chisq=False, note=None):
    
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, is_chisq=is_chisq, note=note)
    s = s + ".predicted-orfs.dna.fa"
    return s

def get_riboseq_predicted_orf_mackowiak_overlap_image(riboseq_base, name, image_type='eps', 
        length=None, offset=None, is_unique=False, is_cds_only=False, is_transcriptome=False, is_chisq=False, note=None):
    
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, is_chisq=is_chisq, note=note)
    s = s + ".predicted-orf-mackowiak-overlap.{}".format(image_type)
    return s

def get_riboseq_predicted_orfs_protein(riboseq_base, name, length=None, offset=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, is_chisq=False, note=None):
    
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, is_chisq=is_chisq, note=note)
    s = s + ".predicted-orfs.protein.fa"
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



def get_riboseq_prob(riboseq_base, name, length=None, is_unique=False, is_cds_only=False, is_transcriptome=False):
    s = get_riboseq_bitseq(riboseq_base, name, length=length, is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome)
    s = s + ".prob"
    return s

def get_riboseq_profiles(riboseq_base, name, length=None, offset=None, is_unique=False, is_cds_only=False, is_transcriptome=False):
    s = get_riboseq_base(riboseq_base, name, 'mtx', length=length, offset=offset, is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome)
    s = s + ".profiles.mtx"
    return s

# r

def get_riboseq_read_filtering_counts(riboseq_base):
    s = os.path.join(riboseq_base, 'read-filtering-counts.csv.gz')
    return s

def get_riboseq_read_filtering_counts_image(riboseq_base, note="", image_type="eps"):
    note_str = ""
    if len(note) > 0:
        note_str = "{}.".format(note)
    fn = "read-filtering-counts.{}{}".format(note_str, image_type)
    s = os.path.join(riboseq_base, fn)
    return s

    


# s
def get_riboseq_sam_base(riboseq_data, name):
    return os.path.join(riboseq_data, 'without-rrna-mapping', 'sam', name)

def get_riboseq_sam_path(riboseq_data):
    return os.path.join(riboseq_data, 'without-rrna-mapping', 'sam')

# t

def get_riboseq_transcript_fasta_file(riboseq_base, name, length=None, is_unique=False, is_cds_only=False, is_transcriptome=False):
    s = get_riboseq_bitseq(riboseq_base, name, length=length, is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome)
    s= s + ".transcripts.fa"
    return s

def get_riboseq_transcript_signal_csv(riboseq_base, name, length=None, is_unique=False, is_cds_only=False, is_transcriptome=False):
    s = get_riboseq_bam_base(riboseq_base, name, length=length, is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome)
    s = s + ".transcript-signal-sums.csv.gz"
    return s


### rna

# b
def get_rna_bam_base(rnaseq_base, name, length=None, is_unique=False, is_cds_only=False, is_transcriptome=False):
    cds_only = get_cds_only_string(is_cds_only)
    unique = get_unique_string(is_unique)
    l = get_length_string(length)
    transcriptome = get_transcriptome_string(is_transcriptome)
    return os.path.join(rnaseq_base, 'mapping', 'bam', '{}{}{}{}{}'.format(name, transcriptome, unique, cds_only, l))


def get_rna_bam(rnaseq_base, name, length=None, is_unique=False, is_cds_only=False, is_transcriptome=False):
    
    s = get_rna_bam_base(rnaseq_base, name, length=length, is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome)
    s = s + ".bam"
    return s

def get_rna_bam_path(base_path):
    return os.path.join(base_path, 'mapping', 'bam')

def get_rna_bam_fastqc_path(base_path):
    return os.path.join(base_path, 'mapping', 'bam', 'fastqc')

def get_rna_bam_fastqc_data(base_path, name):
    fastqc_folder = '{}_name'.format(name)
    return os.path.join(base_path, 'mapping', 'bam', 'fastqc', fastqc_folder, 'fastqc_data.txt')

def get_rna_bitseq(rna_base, name, length=None, is_unique=False, is_cds_only=False, is_transcriptome=False):
    unique = get_unique_string(is_unique)
    cds_only = get_cds_only_string(is_cds_only)
    l = get_length_string(length)
    transcriptome = get_transcriptome_string(is_transcriptome)

    return os.path.join(rna_base, 'mapping', 'transcript-abundance', 'bitseq', '{}-rna{}{}{}{}.bitseq'.format(name, transcriptome, unique, cds_only, l))

def get_rna_bitseq_malphas(rna_base, name, length=None, is_unique=False, is_cds_only=False, is_transcriptome=False):
    s = get_rna_bitseq(rna_base, name, length=length, is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome)
    s = s + ".m_alphas"
    return s
    
# f    
def get_rna_fastq(base_path, name):
    return os.path.join(base_path, 'raw-data', '{}-rna.fastq.gz'.format(name))

# m
def get_rna_most_expressed_isoforms(rna_base, name, length=None, is_unique=False, is_cds_only=False, is_transcriptome=False):
    s = get_rna_bitseq(rna_base, name, length=length, is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome)
    s = s + ".most-expressed-isoforms.csv.gz"
    return s

def get_rna_most_expressed_orfs(rna_base, name, length=None, is_unique=False, is_cds_only=False, is_transcriptome=False):
    s = get_rna_bitseq(rna_base, name, length=length, is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome)
    s = s + ".most-expressed-orfs.csv.gz"
    return s

# p
def get_rna_prob(rna_base, name, length=None, is_unique=False, is_cds_only=False, is_transcriptome=False):
    s = get_rna_bitseq(rna_base, name, length=length, is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome)
    s = s + ".prob"
    return s

# s
def get_rna_sam_base(base_path, name):
    return os.path.join(base_path, 'mapping', 'sam', name)

def get_rna_sam_path(base_path):
    return os.path.join(base_path, 'mapping', 'sam')

# t
def get_rna_transcript_fasta_file(riboseq_base, name, length=None, is_unique=False, is_cds_only=False, is_transcriptome=False):
    s = get_rna_bitseq(riboseq_base, name, length=length, is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome)
    s= s + ".transcripts.fa"
    return s

### s
def get_sequence_features(base_path, name, note=None):
    note_str = ""
    if (note is not None) and  (len(note) > 0):
        note_str = "{}.".format(note)

    fn = '{}.transcript-sequence-features.{}mtx'.format(name, note_str)
    return os.path.join(base_path, 'transcript-index', fn)

# simulation paths

# b
def get_simulation_base(base_path, simulation_index, name, is_micro=False):
    micro = get_micro_string(is_micro)
    s = os.path.join(base_path, 'micropeptide-simulation-data', 'simulation-{}.{}{}'.format(simulation_index, name, micro))
    return s

def get_simulation_bayes_factors(base_path, simulation_index, name):
    s = get_simulation_base(base_path, simulation_index, name)
    s = s + ".bayes-factors.csv.gz"
    return s

def get_simulation_evaluation(base_path, simulation_index, name, is_micro=False):
    s = get_simulation_base(base_path, simulation_index, name, is_micro)
    s = s + ".evaluation.csv.gz"
    return s

# g
def get_simulation_ground_truth(base_path, simulation_index, name):
    s = get_simulation_base(base_path, simulation_index, name)
    s = s + ".ground-truth.csv.gz"
    return s

# o
def get_simulation_orfs(base_path, simulation_index, name):
    s = get_simulation_base(base_path, simulation_index, name)
    s = s + ".transcript-orfs.csv.gz"
    return s

# r
def get_simulation_roc(base_path, simulation_index, name, is_micro=False):
    s = get_simulation_base(base_path, simulation_index, name, is_micro)
    s = s + ".roc.eps"
    return s


# s
def get_simulation_sequence_features(base_path, simulation_index, name):
    s = get_simulation_base(base_path, simulation_index, name)
    s = s + ".transcript-sequence-features.mtx"
    return s

def get_simulation_signals(base_path, simulation_index, name):
    s = get_simulation_base(base_path, simulation_index, name)
    s = s + ".signals.mtx"
    return s

### t
def get_transcript_fasta(base_path, name):
    fn = '{}.transcripts.fa'.format(name)
    return os.path.join(base_path, 'transcript-index', fn)

### u
def get_unannotated_orfs(base_path, name):
    fn = '{}.unannotated-orfs.csv.gz'.format(name)
    return os.path.join(base_path, 'transcript-index', fn)


### w
def get_without_adapters(base_path):
    return os.path.join(base_path, 'without-adapters')

def get_without_adapters_base(base_path, name):
    return os.path.join(base_path, 'without-adapters', name)

def get_without_adapters_fastq(base_path, name):
    return os.path.join(base_path, 'without-adapters', '{}.fastq.gz'.format(name))

def get_without_adapters_fastqc(base_path):
    return os.path.join(base_path, 'without-adapters', 'fastqc')

def get_without_adapters_fastqc_data(base_path, name):
    fastqc_folder = '{}_fastqc'.format(name)
    return os.path.join(base_path, 'without-adapters', 'fastqc', fastqc_folder, 'fastqc_data.txt')


def get_with_rrna(base_path):
    return os.path.join(base_path, 'with-rrna')

def get_with_rrna_fastq(base_path, name):
    return os.path.join(base_path, 'with-rrna', '{}.fastq.gz'.format(name))

def get_with_rrna_fastqc(base_path):
    return os.path.join(base_path, 'with-rrna', 'fastqc')

def get_with_rrna_fastqc_data(base_path, name):
    fastqc_folder = '{}_fastqc'.format(name)
    return os.path.join(base_path, 'with-rrna', 'fastqc', fastqc_folder, 'fastqc_data.txt')

def get_without_rrna(base_path):
    return os.path.join(base_path, 'without-rrna')

def get_without_rrna_fastq(base_path, name):
    return os.path.join(base_path, 'without-rrna', '{}.fastq.gz'.format(name))

def get_without_rrna_fastqc(base_path):
    return os.path.join(base_path, 'without-rrna', 'fastqc')

def get_without_rrna_fastqc_data(base_path, name):
    fastqc_folder = '{}_fastqc'.format(name)
    return os.path.join(base_path, 'without-rrna', 'fastqc', fastqc_folder, 'fastqc_data.txt')


### script descriptions
run_riboseq_preprocessing_description = "This script runs the pipeline developed by Jerry (jxiong@age.mpg.de) for cleaning up ribosome sequencing data. It first runs flexbar to remove adapter sequences added for sequencing. Next, bowtie2 is used to remove reads which map to rRNA. STAR is then used to map the non-rRNA reads to the transcriptome. Finally, samtools is used to create a sorted, indexed bam file. The compression format of the input will be guessed based on its extension (gz or bz2)."

run_rnaseq_preprocessing_description = "This script runs the pipeline developed by Jerry (jxiong@age.mpg.de) for cleaning up totalRNA sequencing data. Flexbar is used to remove adapters and low-quality reads. STAR is used to map the RNA reads to the transcriptome. Finally, samtools is used to create a sorted, indexed bam file. The compression format of the input will be guessed based on its extension (gz or bz2)."

