
"""
    confest.py
"""

import pytest


#REF_DATA_URL = 'https://data.dieterichlab.org/s/JB4AkN7jC574fzp/download'
REF_DATA_URL = 'https://data.dieterichlab.org/s/kYn5sY7YrJPWDiG/download'
REF_LOC = 'c-elegans-chrI-example'
REF_CONFIG = 'c-elegans-test.yaml'


@pytest.fixture(scope="session")
def data_loc(tmp_path_factory):
    
    import pbio.misc.shell_utils as shell_utils
    
    loc = tmp_path_factory.mktemp("data")

    cmd = f"wget --no-verbose {REF_DATA_URL} -O {loc}/data.zip"
    shell_utils.check_call(cmd, call=True, raise_on_error=True)
    
    cmd = f"unzip {loc}/data.zip -d {loc}"
    shell_utils.check_call(cmd, call=True, raise_on_error=True)
    
    return loc


@pytest.fixture(scope="session")
def get_config(data_loc):
    
    from pathlib import Path
    
    loc = Path(data_loc, REF_LOC)
    config = Path(loc, REF_CONFIG)
    config.write_text(config.read_text().replace('/path/to/your/c-elegans-example', loc.as_posix()))
    
    # reference (known) paths from example dataset
    ref_config = {
        'genome_name': 'WBcel235.79.chrI',
        'genome_base_path': Path(loc, 'reference', 'WBcel235.79.chrI').as_posix()
    }

    return (config, ref_config)


@pytest.fixture(scope="session")
def get_genome(get_config):
    
    import pbio.misc.shell_utils as shell_utils
    
    config, ref_config = get_config

    # we could control the output to stdout from running the cmd...
    
    # could we actually test slurm, or other options?
    
    #--num-cpus 6 --mem 50G --use-slurm
    #--logging-level DEBUG --log-file ${LOG}
    cmd = f"prepare-rpbp-genome {config.as_posix()}"
    shell_utils.check_call(cmd, call=True, raise_on_error=True)

    return config, ref_config


@pytest.fixture(scope="session")
def get_files(get_genome):

    import yaml
    import pbio.ribo.ribo_filenames as filenames
    
    config, ref_config = get_genome
    config = yaml.load(open(config), Loader=yaml.FullLoader)
    
    files = {
        'transcript': filenames.get_bed(config['genome_base_path'],
                                        config['genome_name'],
                                        is_annotated=True),
        'fasta': filenames.get_transcript_fasta(config['genome_base_path'],
                                                config['genome_name'],
                                                is_annotated=True),
        'orfs': filenames.get_orfs(config['genome_base_path'],
                                   config['genome_name'],
                                   note=config.get('orf_note'),
                                   is_annotated=True),
        'exons': filenames.get_exons(config['genome_base_path'],
                                     config['genome_name'],
                                     note=config.get('orf_note'),
                                     is_annotated=True),
        'labels': filenames.get_labels(config['genome_base_path'],
                                       config['genome_name'],
                                       note=config.get('orf_note'),
                                       is_annotated=True)
    }
    
    ref_files = {
        'transcript': filenames.get_bed(ref_config['genome_base_path'],
                                        ref_config['genome_name'],
                                        is_annotated=True),
        'fasta': filenames.get_transcript_fasta(ref_config['genome_base_path'],
                                                ref_config['genome_name'],
                                                is_annotated=True),
        'orfs': filenames.get_orfs(ref_config['genome_base_path'],
                                   ref_config['genome_name'],
                                   is_annotated=True),
        'exons': filenames.get_exons(ref_config['genome_base_path'],
                                     ref_config['genome_name'],
                                     is_annotated=True),
        'labels': filenames.get_labels(ref_config['genome_base_path'],
                                       ref_config['genome_name'],
                                       is_annotated=True)
    }

    ret = []
    for key in files:
        ret.append((files[key], ref_files[key]))
    return ret
