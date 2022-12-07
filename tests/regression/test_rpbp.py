"""
    Regression test using the c-elegans reference dataset.
"""

import logging
import pytest
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


def to_df(filename):
    args = {}
    if filename.endswith(".bed.gz") or filename.endswith(".tab.gz"):
        args = {"sep": "\t"}
    elif filename.endswith(".mtx.gz"):
        args = {"sep": " ", "comment": "%", "header": None}
    return pd.read_csv(filename, **args)


# test output of `prepare-rpbp-genome`
def test_pipeline_part1(getf_genome):
    for file, ref_file in getf_genome:
        msg = f"Comparing {file} and {ref_file}"
        logger.info(msg)
        pd.testing.assert_frame_equal(to_df(file), to_df(ref_file), check_exact=True)


# test output of `run-all-rpbp-instances`
@pytest.mark.depends(on=["test_pipeline_part1"])
def test_pipeline_part2(getf_pipeline):
    
    from pbiotools.utils.fastx_utils import get_fasta_dict
    
    # deterministic output - profiles should match exactly
    profiles, periodics, predictions = getf_pipeline
    
    for file, ref_file in profiles:
        msg = f"[PROFILES] Comparing {file} and {ref_file}"
        logger.info(msg)
        pd.testing.assert_frame_equal(to_df(file), to_df(ref_file), check_exact=True)
    
    # periodic-offsets - compare only the periodic lengths
    for file, ref_file in periodics:
        
        file_df = to_df(file[0])
        file_lengths = file[1]
        
        ref_file_df = to_df(ref_file[0])
        ref_file_lengths = ref_file[1]

        file_df = file_df[file_df.length.isin(file_lengths)]
        ref_file_df = ref_file_df[ref_file_df.length.isin(ref_file_lengths)]
        
        msg = f"[PERIODIC-OFFSETS] Comparing {file[0]} and {ref_file[0]} for {', '.join(file_lengths)} " \
              f"and {', '.join(ref_file_lengths)} resp."
        logger.info(msg)
        pd.testing.assert_frame_equal(file_df, ref_file_df, check_exact=False)
    
    # predictions - compare intersection of predictions on deterministic fields only
    cols = ['seqname', 'start', 'end', 'id', 'score', 'strand',
       'thick_start', 'thick_end', 'color', 'num_exons', 'exon_lengths',
       'exon_genomic_relative_starts', 'orf_num', 'orf_len',
       'x_1_sum', 'x_2_sum', 'x_3_sum', 'profile_sum'
    ]
    for file, ref_file in predictions:
        
        file_df = to_df(file)
        ref_file_df = to_df(ref_file)
        
        file_df.columns = file_df.columns.str.replace('#', '')
        ref_file_df.columns = ref_file_df.columns.str.replace('#', '')

        ids = file_df['id'].unique()
        ids_ref = ref_file_df['id'].unique()
        ids_common = np.intersect1d(ids, ids_ref)
        ids_new = np.setdiff1d(ids, ids_ref)
        ids_missed = np.setdiff1d(ids_ref, ids)
        
        file_df = file_df[file_df["id"].isin(ids_common)]
        file_df.reset_index(inplace=True, drop=True)
        ref_file_df = ref_file_df[ref_file_df["id"].isin(ids_common)]
        ref_file_df.reset_index(inplace=True, drop=True)
        
        msg = f"[ORFS BED] Comparing {file} and {ref_file}"
        logger.info(msg)
        pd.testing.assert_frame_equal(file_df[cols], 
                                      ref_file_df[cols], 
                                      check_exact=True)
        
        # alternatively, if we ignore variance estimates, where discrepancies
        # can be very large, then we could compare the means with 
        # some 'reasonable' tolerance... 
        # see https://github.com/dieterich-lab/rp-bp/issues/117
        
        # compare DNA sequences for selected ORFs 
        # proteins are translated using Bio.Seq.translate - don't test
        file_dna = file.replace('.bed.gz', '.dna.fa')
        ref_file_dna = ref_file.replace('.bed.gz', '.dna.fa')
        
        file_dna_dict = get_fasta_dict(file_dna)
        ref_file_dna_dict = get_fasta_dict(ref_file_dna)
        
        file_dna_common = {key: file_dna_dict[key] for key in ids_common}
        ref_file_dna_common = {key: ref_file_dna_dict[key] for key in ids_common}
        
        msg = f"[ORFS DNA SEQS] Comparing {file_dna} and {ref_file_dna}"
        logger.info(msg)
        pd.testing.assert_frame_equal(pd.DataFrame.from_dict(file_dna_common, orient='index'), 
                                      pd.DataFrame.from_dict(ref_file_dna_common, orient='index'), 
                                      check_exact=True)
        
        if ids_new.size != 0:
            msg = f"Unique records predicted due to differences in statistical model " \
                f"estimates: {', '.join(ids_new)}"
            logger.warning(msg)
        if ids_missed.size != 0:
            msg = f"Reference records not predicted: {', '.join(ids_missed)}"
            logger.warning(msg)
        
        
