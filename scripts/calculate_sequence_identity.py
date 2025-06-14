#!/usr/bin/env python3
"""
Calculate Matrix of Sequence Identity
Usage: python sequence_identity.py input.aln > identity_matrix.tsv
"""

import sys
import os
import numpy as np
from collections import OrderedDict

def parse_clustal(filename):
    """Return Sequence Dictionary and Sequence Order"""
    sequences = OrderedDict()
    current_seq = None
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Skip blank lines and annotations
            if not line or line.startswith("CLUSTAL") or line.startswith("#"):
                continue
                
            # Skip lines of conservation markers
            if line.startswith("*") or line.startswith(".") or line.startswith(":"):
                continue
                
            # Analysis lines of sequences
            parts = line.split()
            if not parts:
                continue
                
            seq_id = parts[0]
            seq_data = ''.join(parts[1:-1]) if parts[-1].isdigit() else ''.join(parts[1:])
            
            if seq_id not in sequences:
                sequences[seq_id] = []
            sequences[seq_id].append(seq_data)
    
    # Merge sequence segments
    for seq_id in sequences:
        sequences[seq_id] = ''.join(sequences[seq_id]).upper()
    
    return sequences

def calculate_identity(seq1, seq2):
    """Calculate identity between sequences"""
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length")
    
    valid_positions = 0
    identical = 0
    
    for a, b in zip(seq1, seq2):
        # Ignore blank positions
        if a == '-' and b == '-':
            continue
            
        valid_positions += 1
        if a == b:
            identical += 1
    
    # Avoid errors from zero division
    return identical / valid_positions * 100 if valid_positions > 0 else 0.0

def main():
    if len(sys.argv) != 2:
        print("Error：Please provide input file\nUsage: python sequence_identity.py input.aln")
        sys.exit(1)
    
    input_file = sys.argv[1]
    if not os.path.exists(input_file):
        print(f"Error：File '{input_file}' does not exist")
        sys.exit(1)
    
    # Analysis alignment file
    try:
        sequences = parse_clustal(input_file)
    except Exception as e:
        print(f"Analysis Error: {str(e)}")
        sys.exit(1)
    
    seq_ids = list(sequences.keys())
    n = len(seq_ids)
    
    # Initialize result matrix
    identity_matrix = np.zeros((n, n))
    
    # Calculate sequence identity
    for i in range(n):
        for j in range(i, n):
            identity = calculate_identity(sequences[seq_ids[i]], sequences[seq_ids[j]])
            identity_matrix[i][j] = identity
            identity_matrix[j][i] = identity
    
    # Print result matrix
    print("\t" + "\t".join(seq_ids))  # header
    for i in range(n):
        row = [f"{identity_matrix[i][j]:.2f}" for j in range(n)]
        print(seq_ids[i] + "\t" + "\t".join(row))

if __name__ == "__main__":
    main()
 
