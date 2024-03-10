# %% [markdown]
# # Outline

# %% [markdown]
#   conda create -n ancient_dna_env python=3.8 biopython pandas matplotlib numpy jupyter ipython scipy seaborn -y
# 

# %% [markdown]
# 1. Creating a Directory: Checks if an 'output' directory exists; if not, it creates one.
# 
# 2. Loading Data: Loads AmtDB metadata and sequence IDs from a specified FASTA file.
# 
# 3. Finding Missing Sequences: Identifies sequences present in AmtDB metadata but missing from the FASTA file.
# 
# 4. Loading a Mitochondrial Dataset: Loads a mitochondrial (mt) dataset's sequence IDs and metadata from specified files.
# 
# 5. Extracting and Saving Sequences: Extracts sequences that match specified IDs from a FASTA file and saves them to a new file.
# 
# 6. Matching and Saving Metadata: Matches metadata for specified IDs and saves it to a new CSV file.
# 
# 7. Main Workflow:
#     * Initializes the process by creating directories and loading initial datasets.
#     * Identifies sequences missing in the AmtDB dataset but available in the mitochondrial dataset.
#     * For sequences found, extracts these sequences and their metadata, saving them to the 'output' directory.
#     * Additionally, identifies sequences available in the mitochondrial dataset but not in the AmtDB, extracting and saving these as well.
#     * Completes the execution by indicating the results are in the 'output' directory.

# %% [markdown]
# Create folders 'data/amtDB' and 'data/mitogenomes_reich' and store there the input data which can be downloaded from:
# * https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data
# * https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V54/V54.1.p1/SHARE/public.dir/index_v54.1.p1_MT.html
# * https://amtdb.org/

# %%
from Bio import SeqIO
from collections import defaultdict
import pandas as pd
import os

def create_directories():
    """
    Creates 'output' directory if it does not exist.
    """
    if not os.path.exists('output'):
        print("Creating 'output' directory...\n")
        os.makedirs('output')
    elif not os.path.isdir('data/mitogenomes_reich'):
        print
        # comment: 
    else:
        print("'output' directory already exists.\n")

def load_data(metadata_file, fasta_file):
    """
    Load AmtDB metadata and sequence IDs from the FASTA file.

    Args:
        metadata_file (str): Path to the AmtDB metadata file (e.g. 'amtdb_metadata.csv')
        fasta_file (str): Path to the AmtDB FASTA file (e.g. 'amtdb_1621-samples_7f_a0pkh.fasta')

    Returns:
        meta_amtDB (DataFrame): AmtDB metadata
        ids_seq_fasta (list): List of sequence IDs from the FASTA file (e.g. ['seq1', 'seq2', ...]
    """
    print(f"Loading AmtDB metadata from '{metadata_file}' and sequence IDs from '{fasta_file}'...")
    meta_amtDB = pd.read_csv(metadata_file, sep=',', header=0)
    ids_seq_fasta = [seq_record.id for seq_record in SeqIO.parse(fasta_file, "fasta")]
    print(f"Loaded AmtDB metadata with {len(meta_amtDB)} records and {len(ids_seq_fasta)} sequences.\n")
    return meta_amtDB, ids_seq_fasta

def find_missing_sequences(meta_amtDB, ids_seq_fasta):
    """
    Identifies sequences present in metadata but missing from the FASTA file.

    Args:
        meta_amtDB (DataFrame): AmtDB metadata
        ids_seq_fasta (list): List of sequence IDs from the FASTA file (e.g. ['seq1', 'seq2', ...]

    Returns:
        missing_ids (list): List of sequence IDs present in metadata but missing from the FASTA file
    """
    print("Identifying sequences present in metadata but missing from the FASTA file...")
    amtDB_ids = set(meta_amtDB['identifier'])
    fasta_ids = set(ids_seq_fasta)
    missing_ids = list(amtDB_ids.difference(fasta_ids))
    print(f"Found {len(missing_ids)} missing sequences.\n")
    return missing_ids

def load_mt_dataset(fasta_file, anno_file):
    """
    Load the Reich mt dataset and metadata.

    Args:
        fasta_file (str): Path to the Reich mt dataset FASTA file (e.g. 'mtdna_reich.fasta')
        anno_file (str): Path to the Reich mt dataset metadata file (e.g. 'v54.1.p1_1240K_public.anno')

    Returns:
        ids_mt_dataset (list): List of sequence IDs from the FASTA file (e.g. ['seq1', 'seq2', ...]
        meta_mt_dataset (DataFrame): Reich mt dataset metadata
    """
    print(f"Loading 'Reich mt dataset' from '{fasta_file}' and metadata from '{anno_file}'...")
    ids_mt_dataset = [seq_record.id for seq_record in SeqIO.parse(fasta_file, "fasta")]
    meta_mt_dataset = pd.read_csv(anno_file, sep='\t', header=0, low_memory=False)
    print(f"Loaded 'Reich mt dataset' with {len(ids_mt_dataset)} sequences and metadata with {len(meta_mt_dataset)} records.\n")
    return ids_mt_dataset, meta_mt_dataset

def extract_and_save_sequences(fasta_file, ids, output_file):
    """
    Extracts sequences matching the specified IDs and saves them to a new FASTA file.

    Args:
        fasta_file (str): Path to the input FASTA file
        ids (list): List of sequence IDs to extract
        output_file (str): Path to the output FASTA file
    """
    print(f"Extracting sequences matching the specified IDs from '{fasta_file}'...")
    sequences = [seq_record for seq_record in SeqIO.parse(fasta_file, "fasta") if seq_record.id in ids]
    SeqIO.write(sequences, output_file, "fasta")
    print(f"Saved {len(sequences)} sequences to '{output_file}'.\n")

def match_and_save_metadata(meta_mt_dataset, ids, output_file, id_column_name):
    """
    Matches metadata for the specified IDs and saves it to a new CSV file.

    Args:
        meta_mt_dataset (DataFrame): Reich mt dataset metadata
        ids (list): List of sequence IDs to match
        output_file (str): Path to the output CSV file
        id_column_name (str): Name of the column containing sequence IDs in the metadata
    """
    print(f"Matching metadata for the specified IDs and saving to '{output_file}'...")
    matched_metadata = meta_mt_dataset[meta_mt_dataset[id_column_name].isin(ids)]
    matched_metadata.to_csv(output_file, sep=',', index=False)
    print(f"Saved metadata for {len(matched_metadata)} sequences to '{output_file}'.\n")

def main():
    """
    Main pipeline function.
    """
    create_directories()
    meta_amtDB, ids_seq_fasta_amtDB = load_data('data/amtDB/amtdb_metadata.csv', 'data/amtDB/amtdb_1621-samples_7f_a0pkh.fasta')
    ids_mt_dataset, meta_mt_dataset = load_mt_dataset('data/mitogenomes_reich/mtdna_reich.fasta', 'data/mitogenomes_reich/v54.1.p1_1240K_public/v54.1.p1_1240K_public.anno')

    missing_ids = find_missing_sequences(meta_amtDB, ids_seq_fasta_amtDB)
    available_in_mt_dataset = set(ids_mt_dataset).intersection(set(missing_ids))

    id_column_name = 'Master ID'
    if available_in_mt_dataset:
        extract_and_save_sequences('data/mitogenomes_reich/mtdna_reich.fasta', available_in_mt_dataset, f'output/sequences_missing_in_AmtDB_{len(available_in_mt_dataset)}.fasta')
        match_and_save_metadata(meta_mt_dataset, available_in_mt_dataset, f'output/metadata_for_sequences_missing_in_AmtDB.csv', id_column_name)

    not_in_amtDB = set(ids_mt_dataset) - set(meta_amtDB['identifier'])
    if not_in_amtDB:
        extract_and_save_sequences('data/mitogenomes_reich/mtdna_reich.fasta', not_in_amtDB, f'output/sequences_not_present_in_AmtDB_{len(not_in_amtDB)}.fasta')
        match_and_save_metadata(meta_mt_dataset, not_in_amtDB, f'output/metadata_for_sequences_not_present_in_AmtDB.csv', id_column_name)
        
    print("Pipeline execution complete. Check the 'output' directory for results.\n")

    
    
if __name__ == "__main__":
    main()


# %%
# load the output fasta file and check the number of sequences
# load the output metadata file and check the number of records

missing_seq = SeqIO.parse('output/sequences_missing_in_AmtDB_404.fasta', 'fasta')
missing_seq_metadata = pd.read_csv('output/metadata_for_sequences_missing_in_AmtDB.csv')

not_present_seq = SeqIO.parse('output/sequences_not_present_in_AmtDB_2996.fasta', 'fasta')
not_present_seq_metadata = pd.read_csv('output/metadata_for_sequences_not_present_in_AmtDB.csv')

print(f"Number of missing sequences: {len(list(missing_seq))}")
print(f"Number of metadata records for missing sequences: {len(missing_seq_metadata)}")

print(f"Number of sequences not present in AmtDB: {len(list(not_present_seq))}")
print(f"Number of metadata records for sequences not present in AmtDB: {len(not_present_seq_metadata)}")


# %%
# MASTER ID used for AADR

# Function to load IDs from a CSV file
def load_ids_from_csv(file_path):
    """
    Load IDs from a CSV metadata file.
    """
    df = pd.read_csv(file_path, sep=',', header=0) 
    return set(df['identifier'])

# Function to load IDs from a anno file
def load_ids_from_anno(file_path):
    """
    Load IDs from a anno metadata file.
    """
    df = pd.read_csv(file_path, sep='\t', header=0, low_memory=False)
    return set(df['Master ID'])

# Function to load IDs from a FASTA file
def load_ids_from_fasta(file_path):
    """
    Load IDs from a FASTA file.
    """
    return [seq_record.id for seq_record in SeqIO.parse(file_path, "fasta")]

# Find missing sequences in FASTA given a set of expected IDs
def find_missing_sequences(expected_ids, fasta_ids):
    """
    Identifies expected IDs that are not present in the FASTA IDs.
    """
    return expected_ids - fasta_ids

def extract_and_save_sequences(fasta_file, ids, output_file):
    """
    Extracts sequences matching specified IDs from a FASTA file. In case of duplicate IDs,
    only the sequence with the longest length is kept. If lengths are equal, the first encountered
    sequence is kept. Each ID will correspond to at most one sequence in the output FASTA, 
    ensuring the number of sequences matches the number of unique IDs provided.
    """
    id_to_sequence = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id in ids:
            if record.id not in id_to_sequence:
                id_to_sequence[record.id] = record
            else:
                # If the ID is already present, check if the current sequence is longer
                existing_record = id_to_sequence[record.id]
                if len(record.seq) > len(existing_record.seq):
                    id_to_sequence[record.id] = record
                # If sequences are of the same length, do nothing (keep the existing)
                # This ensures we only update if we find a longer sequence
    
    # Write the unique sequences to the output file
    SeqIO.write(id_to_sequence.values(), output_file, "fasta")
    

# Match metadata for specified IDs and save to CSV AmtDB
def match_and_save_metadata_amtdb(df, ids, output_file):
    """
    Matches metadata for the specified IDs and saves it to a CSV file.
    """
    matched_df = df[df['identifier'].isin(ids)]
    matched_df.to_csv(output_file, index=False)
    
# Match metadata for specified IDs and save to CSV AADR
def match_and_save_metadata_aadr(df, ids, output_file):
    """
    Matches metadata for the specified IDs and saves it to a CSV file.
    """
    matched_df = df[df['Master ID'].isin(ids)]
    matched_df.to_csv(output_file, index=False)
    

# %%
# Load all IDs from both AmtDB and AADR databases.
amtdb_ids_meta = load_ids_from_csv('data/amtDB/amtdb_metadata.csv')
amtdb_ids_fasta = set(load_ids_from_fasta('data/amtDB/amtdb_1621-samples_7f_a0pkh.fasta'))
amtdb_metadata_df = pd.read_csv('data/amtDB/amtdb_metadata.csv', sep=',', header=0)

aadr_ids_meta = load_ids_from_anno('data/mitogenomes_reich/v54.1.p1_1240K_public/v54.1.p1_1240K_public.anno')
aadr_ids_fasta = set(load_ids_from_fasta('data/mitogenomes_reich/mtdna_reich.fasta'))
aadR_metadata_df = pd.read_csv('data/mitogenomes_reich/v54.1.p1_1240K_public/v54.1.p1_1240K_public.anno', sep='\t', header=0, low_memory=False)


# %%
for x, y in zip([amtdb_ids_meta, amtdb_ids_fasta, aadr_ids_meta, aadr_ids_fasta], ['amtdb_ids_meta', 'amtdb_ids_fasta', 'aadr_ids_meta', 'aadr_ids_fasta']):
    print(f"Number of {y}: {len(x)}")
    
print("")

# print first five elements of each set
for x, y in zip([amtdb_ids_meta, amtdb_ids_fasta, aadr_ids_meta, aadr_ids_fasta], ['amtdb_ids_meta', 'amtdb_ids_fasta', 'aadr_ids_meta', 'aadr_ids_fasta']):
    print(f"First five elements of {y}: {list(x)[:5]}")
    

# %%
ids_of_sequences_missing_internally_in_AmtDB = amtdb_ids_meta - amtdb_ids_fasta
ids_of_sequences_not_present_in_AmtDB = aadr_ids_fasta - amtdb_ids_meta  # sequences in AADR but not in AmtDB

print(f"Number of sequences missing internally in AmtDB: {len(ids_of_sequences_missing_internally_in_AmtDB)}")
print(f"Number of sequences not present in AmtDB but present in AADR: {len(ids_of_sequences_not_present_in_AmtDB)}\n")

# save the sequences and metadata
extract_and_save_sequences('data/mitogenomes_reich/mtdna_reich.fasta', ids_of_sequences_missing_internally_in_AmtDB, 'output/sequences_missing_internally_in_AmtDB_masterID.fasta')
match_and_save_metadata_aadr(aadR_metadata_df, ids_of_sequences_missing_internally_in_AmtDB, 'output/metadata_for_sequences_missing_internally_in_AmtDB_masterID.csv')

extract_and_save_sequences('data/mitogenomes_reich/mtdna_reich.fasta', ids_of_sequences_not_present_in_AmtDB, 'output/sequences_not_present_in_AmtDB_masterID.fasta')
match_and_save_metadata_aadr(aadR_metadata_df, ids_of_sequences_not_present_in_AmtDB, 'output/metadata_for_sequences_not_present_in_AmtDB_masterID.csv')

# load the output fasta file and check the number of sequences
# load the output metadata file and check the number of records

missing_seq = SeqIO.parse('output/sequences_missing_internally_in_AmtDB_masterID.fasta', 'fasta')
missing_seq_metadata = pd.read_csv('output/metadata_for_sequences_missing_internally_in_AmtDB_masterID.csv')

not_present_seq = SeqIO.parse('output/sequences_not_present_in_AmtDB_masterID.fasta', 'fasta')
not_present_seq_metadata = pd.read_csv('output/metadata_for_sequences_not_present_in_AmtDB_masterID.csv')

print(f"Number of missing sequences internally in AmtDB, but found in AADR: {len(list(missing_seq))}")
print(f"Number of metadata records for missing sequences: {len(missing_seq_metadata)}\n")

print(f"Number of sequences not present in AmtDB: {len(list(not_present_seq))}")
print(f"Number of metadata records for sequences not present in AmtDB: {len(not_present_seq_metadata)}")


# %%
print(f'Lenght of aadr_ids_fasta: {len(aadr_ids_fasta)}')
print(f'Lenght of amtdb_ids_meta: {len(amtdb_ids_meta)}')
print(f'Lenght of ids_of_sequences_not_present_in_AmtDB: {len(ids_of_sequences_not_present_in_AmtDB)}')

# %%
ids_unique_amtdb = amtdb_ids_meta - aadr_ids_fasta
print(f'Lenght of ids_unique_amtdb: {len(ids_unique_amtdb)}')

# %%
# prunik aadr_ids_fasta a amtdb_ids_meta
ids_common = aadr_ids_fasta.intersection(amtdb_ids_meta)
print(f'Lenght of ids_common: {len(ids_common)}')

# %%
# join two sets together aadr_ids_fasta and amtdb_ids_meta
ids_all = aadr_ids_fasta.union(amtdb_ids_meta)
print(f'Lenght of ids_all: {len(ids_all)}')



# %%
# difference of two sets
ids_unique_aadr = aadr_ids_fasta - amtdb_ids_meta
print(f'Lenght of ids_unique_aadr: {len(ids_unique_aadr)}')
ids_unique_amtdb = amtdb_ids_meta - aadr_ids_fasta
print(f'Lenght of ids_unique_amtdb: {len(ids_unique_amtdb)}')

# %%
ids_all_together = len(ids_unique_aadr) + len(ids_unique_amtdb) + len(ids_common)
print(f'Lenght of ids_all_together: {ids_all_together}')

# %%
print(f'Lenght of aadr_ids_fasta: {len(aadr_ids_fasta)}')
print(f'Lenght of amtdb_ids_meta: {len(amtdb_ids_meta)}')

# %%
not_present_seq = SeqIO.parse('output/sequences_not_present_in_AmtDB_masterID.fasta', 'fasta')
# store the ids in a list and count them
ids_list_not_present = []   
for i, record in enumerate(not_present_seq):
    ids_list_not_present.append(record.id)
print(f"Number of sequences not present in AmtDB: {len(ids_list_not_present)}")
# unique ids
unique_ids_not_present = set(ids_list_not_present)
print(f"Number of unique sequences not present in AmtDB: {len(unique_ids_not_present)}")

# %%
# count the number of records in the metadata file
print(f'Count of records in ids_list_not_present' , len(ids_list_not_present))
# use counter for ids_list_not_present to count the number of occurences of each id
from collections import Counter
counter = Counter(ids_list_not_present)
counter
# find the ones with number > 1
most_common_ids = counter.most_common()
# most_common_ids

# %%
# save to fasta file sequences with the most common ids
most_common_ids_list = [x[0] for x in most_common_ids if x[1] > 1]
most_common_ids_list

# %%

def find_differences(seq1, seq2):
    """Find differences between two sequences."""
    return [i for i in range(min(len(seq1), len(seq2))) if seq1[i] != seq2[i]]

def find_and_compare_duplicates(fasta_file):
    sequences_by_id = defaultdict(list)
    
    for record in SeqIO.parse(fasta_file, 'fasta'):
        sequences_by_id[record.id].append(str(record.seq))
    
    comparison_results = {}
    differences_results = {}  # Store the positions of differences
    
    for id, sequences in sequences_by_id.items():
        if len(sequences) > 1:
            all_identical = all(seq == sequences[0] for seq in sequences)
            comparison_results[id] = all_identical
            if not all_identical:
                # Find differences between the first sequence and the rest
                differences = [find_differences(sequences[0], seq) for seq in sequences[1:]]
                differences_results[id] = differences
    
    return comparison_results, differences_results

fasta_file = 'output/most_common_ids.fasta'
duplicated_sequences, differences_results = find_and_compare_duplicates(fasta_file)

for id, are_sequences_identical in duplicated_sequences.items():
    print(f'ID {id} has identical sequences: {are_sequences_identical}')
    if not are_sequences_identical:
        print(f'   Differences found at positions: {differences_results[id]}\n')
print("")

print(f'Total duplicate IDs: {len(duplicated_sequences)}')
duplicate_count = sum(1 for identical in duplicated_sequences.values() if identical)
print(f'IDs with identical sequences: {duplicate_count}')
print(f'IDs with non-identical sequences: {len(duplicated_sequences) - duplicate_count}')


# %%
extract_and_save_sequences('data/mitogenomes_reich/mtdna_reich.fasta', most_common_ids_list, 'output/most_common_ids.fasta')

# %%
# load metadata_for_sequences_not_present_in_AmtDB_masterID.csv
not_present_seq_metadata = pd.read_csv('output/metadata_for_sequences_not_present_in_AmtDB_masterID.csv')
# not_present_seq_metadata

# %%
"""
Identifying and analyzing the frequency of unique sequence 
IDs not present in AmtDB, ordered by their counts."""

ids_of_sequences_not_present_in_AmtDB
# from ids_of_sequences_not_present_in_AmtDB show unique IDs and their counts
ids_of_sequences_not_present_in_AmtDB = list(ids_of_sequences_not_present_in_AmtDB)
# and their counts
from collections import Counter
ids_of_sequences_not_present_in_AmtDB_counts = Counter(ids_of_sequences_not_present_in_AmtDB)
ids_of_sequences_not_present_in_AmtDB_counts
# order the dictionary by counts
ids_of_sequences_not_present_in_AmtDB_counts = dict(sorted(ids_of_sequences_not_present_in_AmtDB_counts.items(), key=lambda item: item[1], reverse=True))
# show first 5 elements from the dictionary
ids_of_sequences_not_present_in_AmtDB_counts.popitem()

# %% [markdown]
# The difference between number of found metadata for sequences is because of use of 'Master ID' and than 'Genetic ID' in some master ID's is visible .in_preparation added to 'Master ID'.
# And there are more entries of metadata for one sequence. For internally missing is better to use Master ID as identifier, and for sequences not present is better to use Genetic ID as identifier, for better metadata retrieval.

# %% [markdown]
# sequences_missing_internally_in_AmtDB_masterID.fasta
# 
# metadata_for_sequences_missing_internally_in_AmtDB_masterID.csv
# 
# sequences_not_present_in_AmtDB_geneticID.fasta
# 
# metadata_for_sequences_not_present_in_AmtDB_geneticID.csv
# 

# %% [markdown]
# # EIGENSTRAT
# 

# %% [markdown]
# software/EIG/bin/smartpca.perl -i data/mitogenomes_reich/v54.1.p1_HO_public/v54.1.p1_HO_public.geno -a data/mitogenomes_reich/v54.1.p1_HO_public/v54.1.p1_HO_public.snp -b data/mitogenomes_reich/v54.1.p1_HO_public/v54.1.p1_HO_public.ind -k 10 -o output/v54.1.p1_HO_public.pca -p output/v54.1.p1_HO_public.plot -e output/v54.1.p1_HO_public.eval -l output/v54.1.p1_HO_public.log










so we can continue with step number two, revise, optimze and improove the code , provide suggestion how to adapt and change the code to the project outline and where should I put corrected code and how to call it from the ipynb where the analysis will be done. Here is the code.

from Bio import SeqIO
from collections import defaultdict
import pandas as pd
import os

def create_directories():
    """
    Creates 'output' directory if it does not exist.
    """
    if not os.path.exists('output'):
        print("Creating 'output' directory...\n")
        os.makedirs('output')
    elif not os.path.isdir('data/mitogenomes_reich'):
        print
        # comment: 
    else:
        print("'output' directory already exists.\n")

def load_data(metadata_file, fasta_file):
    """
    Load AmtDB metadata and sequence IDs from the FASTA file.

    Args:
        metadata_file (str): Path to the AmtDB metadata file (e.g. 'amtdb_metadata.csv')
        fasta_file (str): Path to the AmtDB FASTA file (e.g. 'amtdb_1621-samples_7f_a0pkh.fasta')

    Returns:
        meta_amtDB (DataFrame): AmtDB metadata
        ids_seq_fasta (list): List of sequence IDs from the FASTA file (e.g. ['seq1', 'seq2', ...]
    """
    print(f"Loading AmtDB metadata from '{metadata_file}' and sequence IDs from '{fasta_file}'...")
    meta_amtDB = pd.read_csv(metadata_file, sep=',', header=0)
    ids_seq_fasta = [seq_record.id for seq_record in SeqIO.parse(fasta_file, "fasta")]
    print(f"Loaded AmtDB metadata with {len(meta_amtDB)} records and {len(ids_seq_fasta)} sequences.\n")
    return meta_amtDB, ids_seq_fasta

def find_missing_sequences(meta_amtDB, ids_seq_fasta):
    """
    Identifies sequences present in metadata but missing from the FASTA file.

    Args:
        meta_amtDB (DataFrame): AmtDB metadata
        ids_seq_fasta (list): List of sequence IDs from the FASTA file (e.g. ['seq1', 'seq2', ...]

    Returns:
        missing_ids (list): List of sequence IDs present in metadata but missing from the FASTA file
    """
    print("Identifying sequences present in metadata but missing from the FASTA file...")
    amtDB_ids = set(meta_amtDB['identifier'])
    fasta_ids = set(ids_seq_fasta)
    missing_ids = list(amtDB_ids.difference(fasta_ids))
    print(f"Found {len(missing_ids)} missing sequences.\n")
    return missing_ids

def load_mt_dataset(fasta_file, anno_file):
    """
    Load the Reich mt dataset and metadata.

    Args:
        fasta_file (str): Path to the Reich mt dataset FASTA file (e.g. 'mtdna_reich.fasta')
        anno_file (str): Path to the Reich mt dataset metadata file (e.g. 'v54.1.p1_1240K_public.anno')

    Returns:
        ids_mt_dataset (list): List of sequence IDs from the FASTA file (e.g. ['seq1', 'seq2', ...]
        meta_mt_dataset (DataFrame): Reich mt dataset metadata
    """
    print(f"Loading 'Reich mt dataset' from '{fasta_file}' and metadata from '{anno_file}'...")
    ids_mt_dataset = [seq_record.id for seq_record in SeqIO.parse(fasta_file, "fasta")]
    meta_mt_dataset = pd.read_csv(anno_file, sep='\t', header=0, low_memory=False)
    print(f"Loaded 'Reich mt dataset' with {len(ids_mt_dataset)} sequences and metadata with {len(meta_mt_dataset)} records.\n")
    return ids_mt_dataset, meta_mt_dataset

def extract_and_save_sequences(fasta_file, ids, output_file):
    """
    Extracts sequences matching the specified IDs and saves them to a new FASTA file.

    Args:
        fasta_file (str): Path to the input FASTA file
        ids (list): List of sequence IDs to extract
        output_file (str): Path to the output FASTA file
    """
    print(f"Extracting sequences matching the specified IDs from '{fasta_file}'...")
    sequences = [seq_record for seq_record in SeqIO.parse(fasta_file, "fasta") if seq_record.id in ids]
    SeqIO.write(sequences, output_file, "fasta")
    print(f"Saved {len(sequences)} sequences to '{output_file}'.\n")

def match_and_save_metadata(meta_mt_dataset, ids, output_file, id_column_name):
    """
    Matches metadata for the specified IDs and saves it to a new CSV file.

    Args:
        meta_mt_dataset (DataFrame): Reich mt dataset metadata
        ids (list): List of sequence IDs to match
        output_file (str): Path to the output CSV file
        id_column_name (str): Name of the column containing sequence IDs in the metadata
    """
    print(f"Matching metadata for the specified IDs and saving to '{output_file}'...")
    matched_metadata = meta_mt_dataset[meta_mt_dataset[id_column_name].isin(ids)]
    matched_metadata.to_csv(output_file, sep=',', index=False)
    print(f"Saved metadata for {len(matched_metadata)} sequences to '{output_file}'.\n")

def main():
    """
    Main pipeline function.
    """
    create_directories()
    meta_amtDB, ids_seq_fasta_amtDB = load_data('data/amtDB/amtdb_metadata.csv', 'data/amtDB/amtdb_1621-samples_7f_a0pkh.fasta')
    ids_mt_dataset, meta_mt_dataset = load_mt_dataset('data/mitogenomes_reich/mtdna_reich.fasta', 'data/mitogenomes_reich/v54.1.p1_1240K_public/v54.1.p1_1240K_public.anno')

    missing_ids = find_missing_sequences(meta_amtDB, ids_seq_fasta_amtDB)
    available_in_mt_dataset = set(ids_mt_dataset).intersection(set(missing_ids))

    id_column_name = 'Master ID'
    if available_in_mt_dataset:
        extract_and_save_sequences('data/mitogenomes_reich/mtdna_reich.fasta', available_in_mt_dataset, f'output/sequences_missing_in_AmtDB_{len(available_in_mt_dataset)}.fasta')
        match_and_save_metadata(meta_mt_dataset, available_in_mt_dataset, f'output/metadata_for_sequences_missing_in_AmtDB.csv', id_column_name)

    not_in_amtDB = set(ids_mt_dataset) - set(meta_amtDB['identifier'])
    if not_in_amtDB:
        extract_and_save_sequences('data/mitogenomes_reich/mtdna_reich.fasta', not_in_amtDB, f'output/sequences_not_present_in_AmtDB_{len(not_in_amtDB)}.fasta')
        match_and_save_metadata(meta_mt_dataset, not_in_amtDB, f'output/metadata_for_sequences_not_present_in_AmtDB.csv', id_column_name)
        
    print("Pipeline execution complete. Check the 'output' directory for results.\n")

    
    
if __name__ == "__main__":
    main()


# %%
# load the output fasta file and check the number of sequences
# load the output metadata file and check the number of records

missing_seq = SeqIO.parse('output/sequences_missing_in_AmtDB_404.fasta', 'fasta')
missing_seq_metadata = pd.read_csv('output/metadata_for_sequences_missing_in_AmtDB.csv')

not_present_seq = SeqIO.parse('output/sequences_not_present_in_AmtDB_2996.fasta', 'fasta')
not_present_seq_metadata = pd.read_csv('output/metadata_for_sequences_not_present_in_AmtDB.csv')

print(f"Number of missing sequences: {len(list(missing_seq))}")
print(f"Number of metadata records for missing sequences: {len(missing_seq_metadata)}")

print(f"Number of sequences not present in AmtDB: {len(list(not_present_seq))}")
print(f"Number of metadata records for sequences not present in AmtDB: {len(not_present_seq_metadata)}")


# %%
# MASTER ID used for AADR

# Function to load IDs from a CSV file
def load_ids_from_csv(file_path):
    """
    Load IDs from a CSV metadata file.
    """
    df = pd.read_csv(file_path, sep=',', header=0) 
    return set(df['identifier'])

# Function to load IDs from a anno file
def load_ids_from_anno(file_path):
    """
    Load IDs from a anno metadata file.
    """
    df = pd.read_csv(file_path, sep='\t', header=0, low_memory=False)
    return set(df['Master ID'])

# Function to load IDs from a FASTA file
def load_ids_from_fasta(file_path):
    """
    Load IDs from a FASTA file.
    """
    return [seq_record.id for seq_record in SeqIO.parse(file_path, "fasta")]

# Find missing sequences in FASTA given a set of expected IDs
def find_missing_sequences(expected_ids, fasta_ids):
    """
    Identifies expected IDs that are not present in the FASTA IDs.
    """
    return expected_ids - fasta_ids

def extract_and_save_sequences(fasta_file, ids, output_file):
    """
    Extracts sequences matching specified IDs from a FASTA file. In case of duplicate IDs,
    only the sequence with the longest length is kept. If lengths are equal, the first encountered
    sequence is kept. Each ID will correspond to at most one sequence in the output FASTA, 
    ensuring the number of sequences matches the number of unique IDs provided.
    """
    id_to_sequence = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id in ids:
            if record.id not in id_to_sequence:
                id_to_sequence[record.id] = record
            else:
                # If the ID is already present, check if the current sequence is longer
                existing_record = id_to_sequence[record.id]
                if len(record.seq) > len(existing_record.seq):
                    id_to_sequence[record.id] = record
                # If sequences are of the same length, do nothing (keep the existing)
                # This ensures we only update if we find a longer sequence
    
    # Write the unique sequences to the output file
    SeqIO.write(id_to_sequence.values(), output_file, "fasta")
    

# Match metadata for specified IDs and save to CSV AmtDB
def match_and_save_metadata_amtdb(df, ids, output_file):
    """
    Matches metadata for the specified IDs and saves it to a CSV file.
    """
    matched_df = df[df['identifier'].isin(ids)]
    matched_df.to_csv(output_file, index=False)
    
# Match metadata for specified IDs and save to CSV AADR
def match_and_save_metadata_aadr(df, ids, output_file):
    """
    Matches metadata for the specified IDs and saves it to a CSV file.
    """
    matched_df = df[df['Master ID'].isin(ids)]
    matched_df.to_csv(output_file, index=False)
    

# %%
# Load all IDs from both AmtDB and AADR databases.
amtdb_ids_meta = load_ids_from_csv('data/amtDB/amtdb_metadata.csv')
amtdb_ids_fasta = set(load_ids_from_fasta('data/amtDB/amtdb_1621-samples_7f_a0pkh.fasta'))
amtdb_metadata_df = pd.read_csv('data/amtDB/amtdb_metadata.csv', sep=',', header=0)

aadr_ids_meta = load_ids_from_anno('data/mitogenomes_reich/v54.1.p1_1240K_public/v54.1.p1_1240K_public.anno')
aadr_ids_fasta = set(load_ids_from_fasta('data/mitogenomes_reich/mtdna_reich.fasta'))
aadR_metadata_df = pd.read_csv('data/mitogenomes_reich/v54.1.p1_1240K_public/v54.1.p1_1240K_public.anno', sep='\t', header=0, low_memory=False)


# %%
for x, y in zip([amtdb_ids_meta, amtdb_ids_fasta, aadr_ids_meta, aadr_ids_fasta], ['amtdb_ids_meta', 'amtdb_ids_fasta', 'aadr_ids_meta', 'aadr_ids_fasta']):
    print(f"Number of {y}: {len(x)}")
    
print("")

# print first five elements of each set
for x, y in zip([amtdb_ids_meta, amtdb_ids_fasta, aadr_ids_meta, aadr_ids_fasta], ['amtdb_ids_meta', 'amtdb_ids_fasta', 'aadr_ids_meta', 'aadr_ids_fasta']):
    print(f"First five elements of {y}: {list(x)[:5]}")
    

# %%
ids_of_sequences_missing_internally_in_AmtDB = amtdb_ids_meta - amtdb_ids_fasta
ids_of_sequences_not_present_in_AmtDB = aadr_ids_fasta - amtdb_ids_meta  # sequences in AADR but not in AmtDB

print(f"Number of sequences missing internally in AmtDB: {len(ids_of_sequences_missing_internally_in_AmtDB)}")
print(f"Number of sequences not present in AmtDB but present in AADR: {len(ids_of_sequences_not_present_in_AmtDB)}\n")

# save the sequences and metadata
extract_and_save_sequences('data/mitogenomes_reich/mtdna_reich.fasta', ids_of_sequences_missing_internally_in_AmtDB, 'output/sequences_missing_internally_in_AmtDB_masterID.fasta')
match_and_save_metadata_aadr(aadR_metadata_df, ids_of_sequences_missing_internally_in_AmtDB, 'output/metadata_for_sequences_missing_internally_in_AmtDB_masterID.csv')

extract_and_save_sequences('data/mitogenomes_reich/mtdna_reich.fasta', ids_of_sequences_not_present_in_AmtDB, 'output/sequences_not_present_in_AmtDB_masterID.fasta')
match_and_save_metadata_aadr(aadR_metadata_df, ids_of_sequences_not_present_in_AmtDB, 'output/metadata_for_sequences_not_present_in_AmtDB_masterID.csv')

# load the output fasta file and check the number of sequences
# load the output metadata file and check the number of records

missing_seq = SeqIO.parse('output/sequences_missing_internally_in_AmtDB_masterID.fasta', 'fasta')
missing_seq_metadata = pd.read_csv('output/metadata_for_sequences_missing_internally_in_AmtDB_masterID.csv')

not_present_seq = SeqIO.parse('output/sequences_not_present_in_AmtDB_masterID.fasta', 'fasta')
not_present_seq_metadata = pd.read_csv('output/metadata_for_sequences_not_present_in_AmtDB_masterID.csv')

print(f"Number of missing sequences internally in AmtDB, but found in AADR: {len(list(missing_seq))}")
print(f"Number of metadata records for missing sequences: {len(missing_seq_metadata)}\n")

print(f"Number of sequences not present in AmtDB: {len(list(not_present_seq))}")
print(f"Number of metadata records for sequences not present in AmtDB: {len(not_present_seq_metadata)}")


# %%
print(f'Lenght of aadr_ids_fasta: {len(aadr_ids_fasta)}')
print(f'Lenght of amtdb_ids_meta: {len(amtdb_ids_meta)}')
print(f'Lenght of ids_of_sequences_not_present_in_AmtDB: {len(ids_of_sequences_not_present_in_AmtDB)}')

# %%
ids_unique_amtdb = amtdb_ids_meta - aadr_ids_fasta
print(f'Lenght of ids_unique_amtdb: {len(ids_unique_amtdb)}')

# %%
# prunik aadr_ids_fasta a amtdb_ids_meta
ids_common = aadr_ids_fasta.intersection(amtdb_ids_meta)
print(f'Lenght of ids_common: {len(ids_common)}')

# %%
# join two sets together aadr_ids_fasta and amtdb_ids_meta
ids_all = aadr_ids_fasta.union(amtdb_ids_meta)
print(f'Lenght of ids_all: {len(ids_all)}')



# %%
# difference of two sets
ids_unique_aadr = aadr_ids_fasta - amtdb_ids_meta
print(f'Lenght of ids_unique_aadr: {len(ids_unique_aadr)}')
ids_unique_amtdb = amtdb_ids_meta - aadr_ids_fasta
print(f'Lenght of ids_unique_amtdb: {len(ids_unique_amtdb)}')

# %%
ids_all_together = len(ids_unique_aadr) + len(ids_unique_amtdb) + len(ids_common)
print(f'Lenght of ids_all_together: {ids_all_together}')

# %%
print(f'Lenght of aadr_ids_fasta: {len(aadr_ids_fasta)}')
print(f'Lenght of amtdb_ids_meta: {len(amtdb_ids_meta)}')

# %%
not_present_seq = SeqIO.parse('output/sequences_not_present_in_AmtDB_masterID.fasta', 'fasta')
# store the ids in a list and count them
ids_list_not_present = []   
for i, record in enumerate(not_present_seq):
    ids_list_not_present.append(record.id)
print(f"Number of sequences not present in AmtDB: {len(ids_list_not_present)}")
# unique ids
unique_ids_not_present = set(ids_list_not_present)
print(f"Number of unique sequences not present in AmtDB: {len(unique_ids_not_present)}")

# %%
# count the number of records in the metadata file
print(f'Count of records in ids_list_not_present' , len(ids_list_not_present))
# use counter for ids_list_not_present to count the number of occurences of each id
from collections import Counter
counter = Counter(ids_list_not_present)
counter
# find the ones with number > 1
most_common_ids = counter.most_common()
# most_common_ids

# %%
# save to fasta file sequences with the most common ids
most_common_ids_list = [x[0] for x in most_common_ids if x[1] > 1]
most_common_ids_list

# %%

def find_differences(seq1, seq2):
    """Find differences between two sequences."""
    return [i for i in range(min(len(seq1), len(seq2))) if seq1[i] != seq2[i]]

def find_and_compare_duplicates(fasta_file):
    sequences_by_id = defaultdict(list)
    
    for record in SeqIO.parse(fasta_file, 'fasta'):
        sequences_by_id[record.id].append(str(record.seq))
    
    comparison_results = {}
    differences_results = {}  # Store the positions of differences
    
    for id, sequences in sequences_by_id.items():
        if len(sequences) > 1:
            all_identical = all(seq == sequences[0] for seq in sequences)
            comparison_results[id] = all_identical
            if not all_identical:
                # Find differences between the first sequence and the rest
                differences = [find_differences(sequences[0], seq) for seq in sequences[1:]]
                differences_results[id] = differences
    
    return comparison_results, differences_results

fasta_file = 'output/most_common_ids.fasta'
duplicated_sequences, differences_results = find_and_compare_duplicates(fasta_file)

for id, are_sequences_identical in duplicated_sequences.items():
    print(f'ID {id} has identical sequences: {are_sequences_identical}')
    if not are_sequences_identical:
        print(f'   Differences found at positions: {differences_results[id]}\n')
print("")

print(f'Total duplicate IDs: {len(duplicated_sequences)}')
duplicate_count = sum(1 for identical in duplicated_sequences.values() if identical)
print(f'IDs with identical sequences: {duplicate_count}')
print(f'IDs with non-identical sequences: {len(duplicated_sequences) - duplicate_count}')


# %%
extract_and_save_sequences('data/mitogenomes_reich/mtdna_reich.fasta', most_common_ids_list, 'output/most_common_ids.fasta')

# %%
# load metadata_for_sequences_not_present_in_AmtDB_masterID.csv
not_present_seq_metadata = pd.read_csv('output/metadata_for_sequences_not_present_in_AmtDB_masterID.csv')
# not_present_seq_metadata

# %%
"""
Identifying and analyzing the frequency of unique sequence 
IDs not present in AmtDB, ordered by their counts."""

ids_of_sequences_not_present_in_AmtDB
# from ids_of_sequences_not_present_in_AmtDB show unique IDs and their counts
ids_of_sequences_not_present_in_AmtDB = list(ids_of_sequences_not_present_in_AmtDB)
# and their counts
from collections import Counter
ids_of_sequences_not_present_in_AmtDB_counts = Counter(ids_of_sequences_not_present_in_AmtDB)
ids_of_sequences_not_present_in_AmtDB_counts
# order the dictionary by counts
ids_of_sequences_not_present_in_AmtDB_counts = dict(sorted(ids_of_sequences_not_present_in_AmtDB_counts.items(), key=lambda item: item[1], reverse=True))
# show first 5 elements from the dictionary
ids_of_sequences_not_present_in_AmtDB_counts.popitem()