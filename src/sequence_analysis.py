#!/usr/bin/env python3
# Author: Jozef Fulop
# Institution: UCT in Prague

from collections import defaultdict
from Bio import SeqIO

class SequenceAnalysis:
    """
    Provides tools for analyzing sequences, including finding differences and identifying duplicates.
    """

    @staticmethod
    def find_differences(seq1, seq2):
        """
        Find differences between two sequences.
        
        Parameters:
            seq1 (str): The first sequence.
            seq2 (str): The second sequence.
            
        Returns:
            list: A list of positions where the sequences differ.
        """
        return [i for i in range(min(len(seq1), len(seq2))) if seq1[i] != seq2[i]]

    def find_and_compare_duplicates(self, fasta_file):
        """
        Identify duplicate sequences in a FASTA file and compare them for differences.
        
        Parameters:
            fasta_file (str): Path to the FASTA file containing sequences.
            
        Returns:
            dict: A dictionary mapping each sequence ID to a boolean indicating whether all its occurrences are identical.
            dict: A dictionary mapping each sequence ID to lists of differences with other occurrences of the same ID.
        """
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
                    differences = [self.find_differences(sequences[0], seq) for seq in sequences[1:]]
                    differences_results[id] = differences

        return comparison_results, differences_results
