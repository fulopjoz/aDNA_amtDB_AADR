#!/usr/bin/env python3
# Author: Jozef Fulop
# Institution: UCT in Prague

from Bio import SeqIO
import pandas as pd


def load_data(metadata_file, fasta_file):
    """
    Load AmtDB metadata and sequence IDs from the FASTA file.

    Args:
        metadata_file (str): Path to the AmtDB metadata file (e.g., 'amtdb_metadata.csv')
        fasta_file (str): Path to the AmtDB FASTA file (e.g., 'amtdb_1621-samples_7f_a0pkh.fasta')

    Returns:
        meta_amtDB (DataFrame): AmtDB metadata.
        ids_seq_fasta (list): List of sequence IDs from the FASTA file.
    """
    print(f"Loading AmtDB metadata from '"
          f"{metadata_file}' and sequence IDs from '{fasta_file}'...")
    
    meta_amtDB = pd.read_csv(metadata_file, sep=',', header=0)
    
    ids_seq_fasta = [
        seq_record.id for seq_record in SeqIO.parse(fasta_file, "fasta")]
    
    print(f"Loaded AmtDB metadata with {len(meta_amtDB)} records and {
          len(ids_seq_fasta)} sequences.\n")
    return meta_amtDB, ids_seq_fasta


def load_mt_dataset(fasta_file, anno_file):
    """
    Load the Reich mt dataset and metadata.

    Args:
        fasta_file (str): Path to the Reich mt dataset FASTA file (e.g., 'mtdna_reich.fasta')
        anno_file (str): Path to the Reich mt dataset metadata file (e.g., 'v54.1.p1_1240K_public.anno')

    Returns:
        ids_mt_dataset (list): List of sequence IDs from the FASTA file.
        meta_mt_dataset (DataFrame): Reich mt dataset metadata.
    """
    print(f"Loading 'Reich mt dataset' from '"
          F"{fasta_file}' and metadata from '{anno_file}'...")
    ids_mt_dataset = [
        seq_record.id for seq_record in SeqIO.parse(fasta_file, "fasta")]
    meta_mt_dataset = pd.read_csv(
        anno_file, sep='\t', header=0, low_memory=False)
    print(f"Loaded 'Reich mt dataset' with {len(
        ids_mt_dataset)} sequences and metadata with {len(meta_mt_dataset)} records.\n")
    return ids_mt_dataset, meta_mt_dataset
