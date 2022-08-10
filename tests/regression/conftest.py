"""
    confest.py
"""

import pytest


# REF_DATA_URL = 'https://data.dieterichlab.org/s/JB4AkN7jC574fzp/download'
REF_DATA_URL = "https://data.dieterichlab.org/s/kYn5sY7YrJPWDiG/download"
REF_LOC = "c-elegans-chrI-example"
REF_CONFIG = "c-elegans-test.yaml"


@pytest.fixture(scope="session")
def data_loc(tmp_path_factory):

    """\
    Download reference dataset for regression testing.

    Parameters
    ----------
    tmp_path_factory
        tmp_path_factory fixture

    Returns
    -------
    :pathlib.Path: base temporary directory
    """

    import pbio.misc.shell_utils as shell_utils

    loc = tmp_path_factory.mktemp("data")

    cmd = f"wget --no-verbose {REF_DATA_URL} -O {loc}/data.zip"
    shell_utils.check_call(cmd, call=True, raise_on_error=True)

    cmd = f"unzip {loc}/data.zip -d {loc}"
    shell_utils.check_call(cmd, call=True, raise_on_error=True)

    return loc


@pytest.fixture(scope="session")
def getf_config(data_loc):

    """\
    Set configuration file.

    Parameters
    ----------
    data_loc
        data_loc fixture

    Returns
    -------
    :obj:`tuple`: configuration files
    """

    import yaml
    from pathlib import Path

    loc = Path(data_loc, REF_LOC)
    config = Path(loc, REF_CONFIG)
    config.write_text(
        config.read_text().replace("/path/to/your/c-elegans-example", loc.as_posix())
    )

    # reference (known) paths from example dataset, keep default options
    ref_config = yaml.load(open(config), Loader=yaml.FullLoader).copy()
    ref_config["genome_name"] = "WBcel235.79.chrI"
    ref_config["genome_base_path"] = Path(
        loc, "reference", "WBcel235.79.chrI"
    ).as_posix()
    ref_config["riboseq_data"] = Path(loc, "reference").as_posix()

    return (config, ref_config)


@pytest.fixture(scope="session")
def get_genome(getf_config):

    """\
    Run `prepare-rpbp-genome`.

    Parameters
    ----------
    getf_config
        getf_config fixture

    Returns
    -------
    :obj:`tuple`: configuration files
    """

    import pbio.misc.shell_utils as shell_utils

    config, ref_config = getf_config

    cmd = f"prepare-rpbp-genome {config.as_posix()}"
    shell_utils.check_call(cmd, call=True, raise_on_error=True)

    return config, ref_config


@pytest.fixture(scope="session")
def getf_genome(get_genome):

    """\
    Get all the Rp-Bp outpout file names
    for the reference genome indices, for the current output
    and the reference dataset.

    Parameters
    ----------
    get_genome
        Fixture calling the index creation

    Returns
    -------
    :obj:`list`: tuples of output files
    """

    import yaml
    import pbio.ribo.ribo_filenames as filenames

    config, ref_config = get_genome
    config = yaml.load(open(config), Loader=yaml.FullLoader)

    files = {
        "transcript": filenames.get_bed(
            config["genome_base_path"], config["genome_name"], is_annotated=True
        ),
        "fasta": filenames.get_transcript_fasta(
            config["genome_base_path"], config["genome_name"], is_annotated=True
        ),
        "orfs": filenames.get_orfs(
            config["genome_base_path"],
            config["genome_name"],
            note=config.get("orf_note"),
            is_annotated=True,
        ),
        "exons": filenames.get_exons(
            config["genome_base_path"],
            config["genome_name"],
            note=config.get("orf_note"),
            is_annotated=True,
        ),
        "labels": filenames.get_labels(
            config["genome_base_path"],
            config["genome_name"],
            note=config.get("orf_note"),
            is_annotated=True,
        ),
    }

    ref_files = {
        "transcript": filenames.get_bed(
            ref_config["genome_base_path"], ref_config["genome_name"], is_annotated=True
        ),
        "fasta": filenames.get_transcript_fasta(
            ref_config["genome_base_path"], ref_config["genome_name"], is_annotated=True
        ),
        "orfs": filenames.get_orfs(
            ref_config["genome_base_path"], ref_config["genome_name"], is_annotated=True
        ),
        "exons": filenames.get_exons(
            ref_config["genome_base_path"], ref_config["genome_name"], is_annotated=True
        ),
        "labels": filenames.get_labels(
            ref_config["genome_base_path"], ref_config["genome_name"], is_annotated=True
        ),
    }

    ret = []
    for key in files:
        ret.append((files[key], ref_files[key]))
    return ret


@pytest.fixture(scope="session")
def get_pipeline(getf_config):

    """\
    Run `run-all-rpbp-instances`.

    Parameters
    ----------
    getf_config
        getf_config fixture

    Returns
    -------
    :obj:`tuple`: configuration files
    """

    import pbio.misc.shell_utils as shell_utils

    config, ref_config = getf_config

    num_cpus = 4
    opts = "--merge-replicates --run-replicates --overwrite --keep-intermediate-files"
    cmd = f"run-all-rpbp-instances {config.as_posix()} " f"--num-cpus {num_cpus} {opts}"
    shell_utils.check_call(cmd, call=True, raise_on_error=True)

    return config, ref_config


@pytest.mark.depends(on=["getf_genome"])
@pytest.fixture(scope="session")
def getf_pipeline(get_pipeline):

    """\
    Get all the Rp-Bp outpout file names
    for the ORF periodicity estimates, ORF profiles,
    and final predictions, for the current output
    and the reference dataset.

    Parameters
    ----------
    get_pipeline
        Fixture calling the pipeline

    Returns
    -------
    :obj:`list`: tuples of output files
    """

    import yaml
    import pbio.ribo.ribo_filenames as filenames
    import pbio.ribo.ribo_utils as ribo_utils

    from rpbp.defaults import metagene_options

    config, ref_config = get_pipeline
    config = yaml.load(open(config), Loader=yaml.FullLoader)

    # identical for config and ref_config
    sample_names = sorted(config["riboseq_samples"].keys())
    riboseq_replicates = ribo_utils.get_riboseq_replicates(config)

    lfiles = [[], []]
    lconfigs = [config, ref_config]

    # we don't test create_base_genome_profile - i.e. output of Flexbar, STAR, etc.
    def populate(s, l, c, merged=False):
        note = c.get("note", None)
        is_unique = not ("keep_riboseq_multimappers" in c)
        fraction = c.get("smoothing_fraction", None)
        reweighting_iterations = c.get("smoothing_reweighting_iterations", None)

        lengths, offsets = None, None
        if not merged:
            lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(
                c,
                s,
                is_unique=is_unique,
                default_params=metagene_options,
            )

            files = {  # part 1 periodicity estimates
                f"metagene_profiles_{s}": filenames.get_metagene_profiles(
                    c["riboseq_data"], s, is_unique=is_unique, note=note
                ),
                f"metagene_profile_bayes_factors_{s}": filenames.get_metagene_profiles_bayes_factors(
                    c["riboseq_data"], s, is_unique=is_unique, note=note
                ),
                f"periodic_offsets_{s}": filenames.get_periodic_offsets(
                    c["riboseq_data"], s, is_unique=is_unique, note=note
                ),
            }
            l.append(files)

        files = {  # part 1 ORF profiles
            f"profiles_{s}": filenames.get_riboseq_profiles(
                c["riboseq_data"],
                s,
                length=lengths,
                offset=offsets,
                is_unique=is_unique,
                note=note,
            ),  # part 2 extended ORF BED (Bayes factor)
            f"bayes_factors_{s}": filenames.get_riboseq_bayes_factors(
                c["riboseq_data"],
                s,
                length=lengths,
                offset=offsets,
                is_unique=is_unique,
                note=note,
                fraction=fraction,
                reweighting_iterations=reweighting_iterations,
            ),
        }
        l.append(files)

        for is_filtered in [True, False]:
            files = {
                f"predicted_orfs_{s}_{is_filtered}": filenames.get_riboseq_predicted_orfs(
                    c["riboseq_data"],
                    s,
                    length=lengths,
                    offset=offsets,
                    is_unique=is_unique,
                    note=note,
                    fraction=fraction,
                    reweighting_iterations=reweighting_iterations,
                    is_filtered=is_filtered,
                ),
                f"predicted_orfs_dna_{s}_{is_filtered}": filenames.get_riboseq_predicted_orfs_dna(
                    c["riboseq_data"],
                    s,
                    length=lengths,
                    offset=offsets,
                    is_unique=is_unique,
                    note=note,
                    fraction=fraction,
                    reweighting_iterations=reweighting_iterations,
                    is_filtered=is_filtered,
                ),
                f"predicted_orfs_protein_{s}_{is_filtered}": filenames.get_riboseq_predicted_orfs_protein(
                    c["riboseq_data"],
                    s,
                    length=lengths,
                    offset=offsets,
                    is_unique=is_unique,
                    note=note,
                    fraction=fraction,
                    reweighting_iterations=reweighting_iterations,
                    is_filtered=is_filtered,
                ),
            }
            l.append(files)

    for s in sample_names:
        for l, c in zip(lfiles, lconfigs):
            populate(s, l, c)

    for s in sorted(riboseq_replicates.keys()):
        for l, c in zip(lfiles, lconfigs):
            populate(s, l, c, merged=True)

    ret = []
    files = {k: v for d in lfiles[0] for k, v in d.items()}
    ref_files = {k: v for d in lfiles[1] for k, v in d.items()}
    for key in files:
        ret.append((files[key], ref_files[key]))
    return ret
