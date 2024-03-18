# tests/test_data_loading.py

import pytest
from src.data_loading import load_data, load_mt_dataset

def test_load_data():
    # Assuming you have a small test dataset for testing purposes
    metadata_file = 'data/amtDB/amtdb_metadata.csv'
    fasta_file = 'data/amtDB/amtdb_1621-samples_7f_a0pkh.fasta'
    
    meta_amtDB, ids_seq_fasta = load_data(metadata_file, fasta_file)
    
    assert not meta_amtDB.empty, "The metadata DataFrame should not be empty."
    assert len(ids_seq_fasta) > 0, "The list of sequence IDs should not be empty."
    assert isinstance(ids_seq_fasta, list), "The sequence IDs should be in a list."

def test_load_mt_dataset():
    # Use the existing dataset for testing
    fasta_file = 'data/mitogenomes_reich/mtdna_reich.fasta'
    anno_file = 'data/mitogenomes_reich/v54.1.p1_1240K_public/v54.1.p1_1240K_public.anno'
    
    ids_mt_dataset, meta_mt_dataset = load_mt_dataset(fasta_file, anno_file)
    
    assert not meta_mt_dataset.empty, "The metadata DataFrame should not be empty."
    assert len(ids_mt_dataset) > 0, "The list of sequence IDs should not be empty."
    assert isinstance(ids_mt_dataset, list), "The sequence IDs should be in a list."





