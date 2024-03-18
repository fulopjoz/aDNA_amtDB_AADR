#!/usr/bin/env python3
# Author: Jozef Fulop
# Institution: UCT in Prague

from Bio import SeqIO

class DataProcessor:
    """
    A class to process and analyze sequence data, focusing on finding missing sequences in one dataset
    and extracting them from another dataset if available.
    """

    def __init__(self, meta_df, fasta_path):
        self.meta_df = meta_df
        self.fasta_path = fasta_path
        self.ids_in_fasta = {
            record.id for record in SeqIO.parse(self.fasta_path, "fasta")}

    def find_missing_sequences(self):
        """
        Identifies sequences present in metadata but missing from the FASTA file associated with this processor.
        """
        meta_ids = set(self.meta_df['identifier'])
        missing_ids = list(meta_ids.difference(self.ids_in_fasta))
        print(f"Found {len(missing_ids)} missing sequences.\n")
        return missing_ids

    def extract_and_save_sequences(self, fasta_file, ids, output_file):
        """
        Extracts sequences matching specified IDs from a given FASTA file and saves them to a new FASTA file.

        Args:
            fasta_file (str): Path to the FASTA file from which to extract sequences.
            ids (list): List of sequence IDs to extract.
            output_file (str): Path to the output FASTA file where extracted sequences will be saved.
        """
        print(f"Extracting sequences matching the specified IDs from '{
              fasta_file}'...")
        sequences = [record for record in SeqIO.parse(
            fasta_file, "fasta") if record.id in ids]
        SeqIO.write(sequences, output_file, "fasta")
        print(f"Saved {len(sequences)} sequences to '{output_file}'.\n")

    def validate_ids(self):
        """
        Validates the overlap between sequence IDs in the metadata and the FASTA file associated with this processor.
        """
        meta_ids = set(self.meta_df['identifier'])
        intersection_ids = meta_ids.intersection(self.ids_in_fasta)
        if not intersection_ids:
            print("Warning: No overlapping IDs found between metadata and FASTA file.")
        else:
            print(f"Validation: Found {
                  len(intersection_ids)} overlapping IDs.")


class MetadataMatcher:
    """
    A class designed to match and manage metadata for specific sequence IDs.
    """
    @staticmethod
    def match_and_save_metadata(df, ids, output_file, id_column_name="identifier"):
        """
        Matches metadata for specified IDs and saves it to a CSV file.
        """
        matched_metadata = df[df[id_column_name].isin(ids)]
        matched_metadata.to_csv(output_file, sep=',', index=False)
        print(f"Saved metadata for {
              len(matched_metadata)} sequences to '{output_file}'.\n")
