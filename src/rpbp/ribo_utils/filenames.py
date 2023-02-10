import glob
import os


def get_annotated_string(is_annotated):
    annotated = ""
    if is_annotated:
        annotated = ".annotated"
    return annotated


def get_chisq_string(is_chisq):
    chisq = ""
    if is_chisq:
        chisq = ".chisq"
    return chisq


def get_de_novo_string(is_de_novo):
    de_novo = ""
    if is_de_novo:
        de_novo = ".de-novo"
    return de_novo


def get_filtered_string(is_filtered):
    filtered_str = ""
    if is_filtered:
        filtered_str = ".filtered"
    return filtered_str


def get_fraction_string(fraction=None):
    fraction_str = ""
    if (fraction is not None) and (len(str(fraction)) > 0):
        fraction_str = ".frac-{}".format(fraction)
    return fraction_str


def get_length_string(length=None):
    l = ""
    if length is not None:
        if isinstance(length, (list, tuple)):
            l = "-".join(str(l) for l in length)
            l = ".length-{}".format(l)
        else:
            l = ".length-{}".format(str(length))
    return l


def get_note_string(note=None):
    note_str = ""
    if (note is not None) and (len(note) > 0):
        note_str = ".{}".format(note)
    return note_str


def get_offset_string(offset=None):
    o = ""
    if offset is not None:
        if isinstance(offset, (list, tuple)):
            o = "-".join(str(o) for o in offset)
            o = ".offset-{}".format(o)
        else:
            o = ".offset-{}".format(str(offset))
    return o


def get_reweighting_iterations_string(reweighting_iterations=None):
    reweighting_iterations_str = ""
    if (reweighting_iterations is not None) and (len(str(reweighting_iterations)) > 0):
        reweighting_iterations_str = ".rw-{}".format(reweighting_iterations)
    return reweighting_iterations_str


def get_star_input_string(is_star_input):
    s = ""
    if is_star_input:
        s = ".star-input"
    return s


def get_unique_string(is_unique):
    unique = ""
    if is_unique:
        unique = "-unique"
    return unique


### filenames


# b


def get_riboseq_base(
    riboseq_base,
    name,
    sub_folder,
    length=None,
    offset=None,
    is_unique=False,
    fraction=None,
    reweighting_iterations=None,
    is_chisq=False,
    is_filtered=False,
    note=None,
):

    unique = get_unique_string(is_unique)
    l = get_length_string(length)
    o = get_offset_string(offset)
    chisq = get_chisq_string(is_chisq)
    n = get_note_string(note)
    f = get_fraction_string(fraction)
    r = get_reweighting_iterations_string(reweighting_iterations)
    fi = get_filtered_string(is_filtered)

    fn = "".join([name, n, unique, l, o, f, r, chisq, fi])

    return os.path.join(riboseq_base, sub_folder, fn)


def get_riboseq_bam_base(riboseq_base, name, **kwargs):

    bam_base = get_riboseq_base(riboseq_base, name, "without-rrna-mapping", **kwargs)

    return bam_base


def get_riboseq_bam(riboseq_base, name, **kwargs):

    s = get_riboseq_bam_base(riboseq_base, name, **kwargs)
    s = s + ".bam"
    return s


def get_riboseq_bam_fastqc_path(riboseq_data):
    return os.path.join(riboseq_data, "without-rrna-mapping", "fastqc")


def get_riboseq_bam_fastqc_data(
    riboseq_data,
    name,
    length=None,
    is_unique=False,
    note=None,
    is_chisq=False,
):

    unique = get_unique_string(is_unique)
    l = get_length_string(length)
    n = get_note_string(note)
    c = get_chisq_string(is_chisq)
    name = "{}{}{}{}{}".format(name, n, unique, l, c)

    fastqc_folder = "{}_fastqc".format(name)
    return os.path.join(
        riboseq_data, "without-rrna-mapping", "fastqc", fastqc_folder, "fastqc_data.txt"
    )


def get_bed(
    base_path,
    name,
    is_annotated=False,
    is_de_novo=False,
):

    c = get_annotated_string(is_annotated)
    d = get_de_novo_string(is_de_novo)
    fn = "{}{}{}.bed.gz".format(name, c, d)
    return os.path.join(base_path, fn)


def get_riboseq_bayes_factors(riboseq_base, name, **kwargs):

    s = get_riboseq_base(riboseq_base, name, "orf-predictions", **kwargs)

    s = s + ".bayes-factors.bed.gz"
    return s


# c


def get_riboseq_cell_type_protein(riboseq_base, name, **kwargs):

    s = get_riboseq_base(riboseq_base, name, "cell-types", **kwargs)
    s = s + ".predicted-orfs.protein.fa"
    return s


# e


def get_exons(base_path, name, is_annotated=False, is_de_novo=False, note=None):

    note_str = get_note_string(note)
    c = get_annotated_string(is_annotated)
    d = get_de_novo_string(is_de_novo)
    fn = "{}.orfs-exons{}{}{}.bed.gz".format(name, c, d, note_str)
    return os.path.join(base_path, "transcript-index", fn)


def get_labels(base_path, name, is_annotated=False, is_de_novo=False, note=None):

    note_str = get_note_string(note)
    c = get_annotated_string(is_annotated)
    d = get_de_novo_string(is_de_novo)
    fn = "{}.orfs-labels{}{}{}.tab.gz".format(name, c, d, note_str)
    return os.path.join(base_path, "transcript-index", fn)


# f


def get_fastqc_name(filename):
    """Given the sequence or alignment filename, this function extracts the name
    used by fastqc. In particular, it removes the ".fasta[.gz]" or ".fastq[.gz]"
    or ".sam" or ".bam" ending from the filename.

    Args:
        filename: the name of the sequence or alignment file (NOT including the path)

    Returns:
        string : the fastqc name

    Imports:
        misc.utils
    """

    import pbiotools.misc.utils as utils

    # first, get the filename
    filename = utils.get_basename(filename)
    filename = filename.replace(".fasta", "")
    filename = filename.replace(".fastq", "")
    filename = filename.replace(".sam", "")
    filename = filename.replace(".bam", "")
    filename = filename.replace(".gz", "")

    return filename


def get_riboseq_fastq(riboseq_data, name):
    return os.path.join(riboseq_data, "raw-data", "{}.fastq.gz".format(name))


def get_riboseq_frame_counts(riboseq_base, name, **kwargs):
    sub_folder = kwargs.pop("sub_folder", "")
    base = get_riboseq_base(riboseq_base, name, sub_folder, **kwargs)
    loc = f"{base}.frame-counts.csv.gz"
    return loc


# g


def get_gtf(
    config,
):

    if "de_novo_gtf" in config:
        base_path = config["genome_base_path"]
        name = config["genome_name"]
        fn = f"{name}.gtf"
        return os.path.join(base_path, fn)
    else:
        return config["gtf"]


# m


def get_metagene_profile_image(base, name, image_type="eps", **kwargs):

    s = get_riboseq_base(base, name, "metagene-profiles", **kwargs)
    s = s + "." + image_type
    return s


def get_metagene_profiles(riboseq_base, name, **kwargs):

    s = get_riboseq_base(riboseq_base, name, "metagene-profiles", **kwargs)
    s = s + ".metagene-profile.csv.gz"
    return s


def get_metagene_profiles_bayes_factors(riboseq_base, name, **kwargs):

    s = get_riboseq_base(riboseq_base, name, "metagene-profiles", **kwargs)
    s = s + ".metagene-periodicity-bayes-factors.csv.gz"
    return s


def get_default_models_base():
    import os
    import inspect
    import rpbp

    return os.path.join(os.path.dirname(inspect.getfile(rpbp)), "models")


def get_models(models_base, model_type):
    import shlex

    # query via stan_file
    path_ex = os.path.join(models_base, model_type, "*stan")
    models = glob.glob(path_ex)
    models = [shlex.quote(m) for m in models]
    return models


def get_metagene_profile_bayes_factor_image(
    riboseq_base, name, image_type="eps", **kwargs
):

    s = get_riboseq_base(riboseq_base, name, "metagene-profiles", **kwargs)
    s = s + ".bayes-factors." + image_type
    return s


# o


def get_orfs(base_path, name, is_annotated=False, is_de_novo=False, note=None):
    note_str = get_note_string(note)
    c = get_annotated_string(is_annotated)
    d = get_de_novo_string(is_de_novo)
    fn = "{}.orfs-genomic{}{}{}.bed.gz".format(name, c, d, note_str)
    return os.path.join(base_path, "transcript-index", fn)


def get_orf_type_profile_base(
    riboseq_base,
    name,
    length=None,
    offset=None,
    is_unique=False,
    fraction=None,
    reweighting_iterations=None,
    note=None,
    is_chisq=False,
    subfolder="orf-predictions",
):

    subfolder = os.path.join(subfolder, "plots")
    s = get_riboseq_base(
        riboseq_base,
        name,
        subfolder,
        length=length,
        offset=offset,
        is_unique=is_unique,
        fraction=fraction,
        reweighting_iterations=reweighting_iterations,
        note=note,
        is_chisq=is_chisq,
    )
    return s


def get_orf_type_profile_image(base_path, orf_type, strand, image_type="eps"):
    fn = ".{}.{}.metagene-profiles.{}".format(orf_type, strand, image_type)
    return base_path + fn


# p


def get_periodic_offsets(riboseq_base, name, **kwargs):

    sub_folder = kwargs.pop("sub_folder", "metagene-profiles")
    s = get_riboseq_base(riboseq_base, name, sub_folder, **kwargs)
    s = s + ".periodic-offsets.csv.gz"
    return s


def get_riboseq_peptide_matches(riboseq_base, name, peptide_name, **kwargs):

    n = "{}-{}".format(name, peptide_name)

    s = get_riboseq_base(riboseq_base, n, "peptide-matches", **kwargs)
    s = s + ".peptide-matches.csv.gz"
    return s


def get_riboseq_predicted_orfs(riboseq_base, name, **kwargs):
    sub_folder = kwargs.pop("sub_folder", "orf-predictions")
    base = get_riboseq_base(riboseq_base, name, sub_folder, **kwargs)
    loc = f"{base}.predicted-orfs.bed.gz"
    return loc


def get_riboseq_predicted_orfs_dna(riboseq_base, name, **kwargs):

    s = get_riboseq_base(riboseq_base, name, "orf-predictions", **kwargs)
    s = s + ".predicted-orfs.dna.fa"
    return s


def get_riboseq_predicted_orfs_protein(riboseq_base, name, **kwargs):

    s = get_riboseq_base(riboseq_base, name, "orf-predictions", **kwargs)
    s = s + ".predicted-orfs.protein.fa"
    return s


def get_riboseq_profiles(riboseq_base, name, **kwargs):

    s = get_riboseq_base(riboseq_base, name, "orf-profiles", **kwargs)
    s = s + ".profiles.mtx.gz"
    return s


# r


def get_riboseq_read_filtering_counts(riboseq_base, name, **kwargs):
    sub_folder = kwargs.pop("sub_folder", "")
    base = get_riboseq_base(riboseq_base, name, sub_folder, **kwargs)
    loc = f"{base}.read-filtering-counts.csv.gz"
    return loc


def get_riboseq_read_length_distribution(riboseq_base, name, **kwargs):
    sub_folder = kwargs.pop("sub_folder", "without-rrna-mapping")
    base = get_riboseq_base(riboseq_base, name, sub_folder, **kwargs)
    loc = f"{base}.length-distribution.csv.gz"
    return loc


def get_raw_data_path(base_path):
    return os.path.join(base_path, "raw-data")


def get_raw_data_fastqc_path(base_path):
    rdp = get_raw_data_path(base_path)
    return os.path.join(rdp, "fastqc")


def get_raw_data_fastqc_data(base_path, filename):

    name = get_fastqc_name(filename)
    fastqc_folder = "{}_fastqc".format(name)
    rdp = get_raw_data_fastqc_path(base_path)

    p = os.path.join(rdp, fastqc_folder, "fastqc_data.txt")
    return p


# s


def get_star_index(base_path, name):
    fn = "{}".format(name)
    return os.path.join(base_path, "STAR", fn)


# t


def get_transcript_fasta(
    base_path,
    name,
    is_annotated=False,
    is_de_novo=False,
):

    c = get_annotated_string(is_annotated)
    d = get_de_novo_string(is_de_novo)
    fn = "{}.transcripts{}{}.fa".format(name, c, d)
    return os.path.join(base_path, "transcript-index", fn)


# w


def get_without_adapters_base(base_path, name, note=None):
    n = get_note_string(note)
    base = "{}{}".format(name, n)
    return os.path.join(base_path, "without-adapters", base)


def get_without_adapters_fastq(base_path, name, note=None):
    base = get_without_adapters_base(base_path, name, note=note)
    fastq = "{}.fastq.gz".format(base)
    return fastq


def get_without_adapters_fastqc(base_path):
    return os.path.join(base_path, "without-adapters", "fastqc")


def get_without_adapters_fastqc_data(base_path, name, note=None):
    n = get_note_string(note)
    fastqc_folder = "{}{}_fastqc".format(name, n)
    return os.path.join(
        base_path, "without-adapters", "fastqc", fastqc_folder, "fastqc_data.txt"
    )


def get_with_rrna_fastq(base_path, name, note=None):
    n = get_note_string(note)
    name = "{}{}".format(name, n)
    return os.path.join(base_path, "with-rrna", "{}.fastq.gz".format(name))


def get_with_rrna_fastqc(base_path):
    return os.path.join(base_path, "with-rrna", "fastqc")


def get_with_rrna_fastqc_data(base_path, name, note=None):
    n = get_note_string(note)
    fastqc_folder = "{}{}_fastqc".format(name, n)
    return os.path.join(
        base_path, "with-rrna", "fastqc", fastqc_folder, "fastqc_data.txt"
    )


def get_without_rrna_fastq(base_path, name, note=None):
    n = get_note_string(note)
    name = "{}{}".format(name, n)
    return os.path.join(base_path, "without-rrna", "{}.fastq.gz".format(name))


def get_without_rrna_fastqc(base_path):
    return os.path.join(base_path, "without-rrna", "fastqc")


def get_without_rrna_fastqc_data(base_path, name, note=None):
    n = get_note_string(note)
    fastqc_folder = "{}{}_fastqc".format(name, n)
    return os.path.join(
        base_path, "without-rrna", "fastqc", fastqc_folder, "fastqc_data.txt"
    )
