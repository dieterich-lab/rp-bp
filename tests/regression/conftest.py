"""
confest.py
"""

import pytest
import sys

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
    import pbiotools.misc.shell_utils as shell_utils

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
    import pbiotools.misc.shell_utils as shell_utils

    config, ref_config = getf_config

    num_cpus = 6
    star_options = '--star-options "--genomeSAindexNbases 10"'
    cmd = (
        f"prepare-rpbp-genome {config.as_posix()} "
        f"--num-cpus {num_cpus} {star_options}"
    )
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
    import rpbp.ribo_utils.filenames as filenames

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
        ),
        "exons": filenames.get_exons(
            config["genome_base_path"],
            config["genome_name"],
            note=config.get("orf_note"),
        ),
        "labels": filenames.get_labels(
            config["genome_base_path"],
            config["genome_name"],
            note=config.get("orf_note"),
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
            ref_config["genome_base_path"], ref_config["genome_name"]
        ),
        "exons": filenames.get_exons(
            ref_config["genome_base_path"], ref_config["genome_name"]
        ),
        "labels": filenames.get_labels(
            ref_config["genome_base_path"], ref_config["genome_name"]
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
    import pbiotools.misc.shell_utils as shell_utils

    config, ref_config = getf_config

    if sys.platform == "darwin":
        num_cpus = 1  # avoid parallel processing issue on macos: https://github.com/dieterich-lab/rp-bp/issues/140
    else:
        num_cpus = 6  # multiprocessing.cpu_count() see https://github.com/dieterich-lab/rp-bp/issues/144
    opts = (
        "--merge-replicates --run-replicates --overwrite "
        "--keep-intermediate-files --write-unfiltered"
    )
    cmd = f"run-all-rpbp-instances {config.as_posix()} --num-cpus {num_cpus} {opts}"
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

    Note: Some files are not used for testing, but
    we leave them here to have a track record of all
    output files.

    Parameters
    ----------
    get_pipeline
        Fixture calling the pipeline

    Returns
    -------
    :obj:`list`: tuples of output files
    """

    import yaml
    import rpbp.ribo_utils.utils as ribo_utils
    import rpbp.ribo_utils.filenames as filenames

    from rpbp.defaults import metagene_options

    config, ref_config = get_pipeline
    config = yaml.load(open(config), Loader=yaml.FullLoader)

    # identical for config and ref_config
    sample_names = sorted(config["riboseq_samples"].keys())
    riboseq_replicates = ribo_utils.get_riboseq_replicates(config)

    lfiles = [[], []]
    lconfigs = [config, ref_config]

    # we don't test create_base_genome_profile - i.e. output of Flexbar, STAR, etc.
    def populate(name, lf, lc, merged=False):
        note = lc.get("note", None)
        is_unique = not ("keep_riboseq_multimappers" in lc)
        fraction = lc.get("smoothing_fraction", None)
        reweighting_iterations = lc.get("smoothing_reweighting_iterations", None)

        lengths, offsets = None, None
        if not merged:
            lengths, offsets = ribo_utils.get_periodic_lengths_and_offsets(
                lc,
                name,
                is_unique=is_unique,
                default_params=metagene_options,
            )

            files = {  # part 1 periodicity estimates
                f"metagene_profiles_{name}": filenames.get_metagene_profiles(
                    lc["riboseq_data"], name, is_unique=is_unique, note=note
                ),  # metagene_profile_bayes_factors UNUSED - see periodic_offsets
                f"metagene_profile_bayes_factors_{name}": filenames.get_metagene_profiles_bayes_factors(
                    lc["riboseq_data"], name, is_unique=is_unique, note=note
                ),  # only comparing periodic lengths
                f"periodic_offsets_{name}": (
                    filenames.get_periodic_offsets(
                        lc["riboseq_data"], name, is_unique=is_unique, note=note
                    ),
                    lengths,
                ),
            }
            lf.append(files)

        files = {  # part 1 ORF profiles
            f"profiles_{name}": filenames.get_riboseq_profiles(
                lc["riboseq_data"],
                name,
                length=lengths,
                offset=offsets,
                is_unique=is_unique,
                note=note,
            ),  # part 2 extended ORF BED (bayes_factors UNUSED -  see predicted_orfs)
            f"bayes_factors_{name}": filenames.get_riboseq_bayes_factors(
                lc["riboseq_data"],
                name,
                length=lengths,
                offset=offsets,
                is_unique=is_unique,
                note=note,
                fraction=fraction,
                reweighting_iterations=reweighting_iterations,
            ),
        }
        lf.append(files)

        # test all using [--write-unfiltered]
        for is_filtered in [True, False]:
            files = {
                f"predicted_orfs_{name}_{is_filtered}": filenames.get_riboseq_predicted_orfs(
                    lc["riboseq_data"],
                    name,
                    length=lengths,
                    offset=offsets,
                    is_unique=is_unique,
                    note=note,
                    fraction=fraction,
                    reweighting_iterations=reweighting_iterations,
                    is_filtered=is_filtered,
                ),  # UNUSED from list, but tested - see test_rpbp.py
                f"predicted_orfs_dna_{name}_{is_filtered}": filenames.get_riboseq_predicted_orfs_dna(
                    lc["riboseq_data"],
                    name,
                    length=lengths,
                    offset=offsets,
                    is_unique=is_unique,
                    note=note,
                    fraction=fraction,
                    reweighting_iterations=reweighting_iterations,
                    is_filtered=is_filtered,
                ),  # UNUSED
                f"predicted_orfs_protein_{name}_{is_filtered}": filenames.get_riboseq_predicted_orfs_protein(
                    lc["riboseq_data"],
                    name,
                    length=lengths,
                    offset=offsets,
                    is_unique=is_unique,
                    note=note,
                    fraction=fraction,
                    reweighting_iterations=reweighting_iterations,
                    is_filtered=is_filtered,
                ),
            }
            lf.append(files)

    for name in sample_names:
        for lf, lc in zip(lfiles, lconfigs):
            populate(name, lf, lc)

    for name in sorted(riboseq_replicates.keys()):
        for lf, lc in zip(lfiles, lconfigs):
            populate(name, lf, lc, merged=True)

    ret_exact, ret_periodic, ret_predictions = [], [], []
    files = {k: v for d in lfiles[0] for k, v in d.items()}
    ref_files = {k: v for d in lfiles[1] for k, v in d.items()}

    for key in files:
        if "metagene_profiles" in key or "profiles_" in key:
            ret_exact.append((files[key], ref_files[key]))
        if "periodic_offsets_" in key:
            ret_periodic.append((files[key], ref_files[key]))
        if "predicted_orfs" in key and "dna" not in key and "protein" not in key:
            ret_predictions.append((files[key], ref_files[key]))

    return (ret_exact, ret_periodic, ret_predictions)
