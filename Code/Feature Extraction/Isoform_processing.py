# GOALS 
# (1) To map isoform IDs to protein-coding gene IDs, extracted from FlyBase genotype_phenotype.tsv (FB2025_02)
# (2) To identify the number of isoforms per protein-coding gene in each phenotypic class  
# (3) To extract CDS sequences for the longest coding isoform of each labelled gene 

# Load libraries
import pandas as pd
import numpy as np
import os
from collections import defaultdict
import pandas as pd

# Define function to get the number of isoforms per protein-coding gene
from collections import defaultdict
import pandas as pd

def parse_fasta_with_lengths_from_header(filepath):
    """
    Parses a FASTA file with FlyBase-style CDS headers to:
    - Map each FBgn (gene) ID to its isoform FBpp IDs
    - Use the 'length=' field in the header as the protein sequence length

    Parameters:
    -----------
    filepath : str
        Path to the FASTA file.

    Returns:
    --------
    tuple:
        - dict: FBgn_ID -> list of FBpp_IDs
        - dict: FBpp_ID -> sequence length (from header)
    """
    parent_to_isoforms = defaultdict(list)
    isoform_lengths = {}

    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('>'):
                header = line[1:].strip()
                parts = header.split()

                # Extract FBpp ID
                fbpp_id = parts[0]

                # Extract parent FBgn ID(s)
                parent_field = next((p for p in parts if p.startswith("parent=")), "")
                parents = parent_field.split('=')[1].split(',') if parent_field else []
                fbgn_parents = [p for p in parents if p.startswith("FBgn")]

                # Extract length=
                length_field = next((p for p in parts if p.startswith("length=")), "")
                length_val_str = length_field.split('=')[1].rstrip(';')
                length_val = int(length_val_str) if length_val_str.isdigit() else 0

                # Map FBgn → FBpp
                for fbgn in fbgn_parents:
                    parent_to_isoforms[fbgn].append(fbpp_id)

                # Store length
                isoform_lengths[fbpp_id] = length_val

    return parent_to_isoforms, isoform_lengths

def parse_fasta_with_lengths_from_header(filepath):
    parent_to_isoforms = defaultdict(list)
    isoform_lengths = {}

    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('>'):
                header = line[1:].strip()
                parts = header.split()

                # Extract FBpp ID
                fbpp_id = parts[0]

                # Initialize empty defaults
                fbgn_parents = []
                length_val = 0

                for p in parts[1:]:
                    if p.startswith("parent="):
                        parent_str = p.split('=', 1)[1].strip(';')
                        parents = parent_str.split(',')
                        fbgn_parents = [x for x in parents if x.startswith("FBgn")]
                    elif p.startswith("length="):
                        length_str = p.split('=', 1)[1].strip(';')
                        if length_str.isdigit():
                            length_val = int(length_str)

                # Map FBgn → FBpp
                for fbgn in fbgn_parents:
                    parent_to_isoforms[fbgn].append(fbpp_id)

                # Map FBpp → length
                isoform_lengths[fbpp_id] = length_val

    return parent_to_isoforms, isoform_lengths

# === Usage ===

fasta_file = '/Users/tanishi/FYP/CDS/NEG_CDS_4.fasta'

parent_to_isoforms, isoform_lengths = parse_fasta_with_lengths_from_header(fasta_file)

# Build DataFrame
data = []
for fbgn, isoforms in parent_to_isoforms.items():
    # Find the longest isoform based on length
    longest_isoform = max(
        isoforms,
        key=lambda fbpp: isoform_lengths.get(fbpp, 0)
    )
    longest_isoform_length = isoform_lengths.get(longest_isoform, 0)
    
    data.append({
        'FBgn_ID': fbgn,
        'Num_Isoforms': len(isoforms),
        'Longest_Isoform_ID': longest_isoform,
        'Longest_Isoform_Length': longest_isoform_length
    })

NEG_4_isoforms = pd.DataFrame(data)
NEG_4_isoforms['FBgn_ID'].nunique()
print(NEG_4_isoforms.head())

NEG_4_isoforms.to_csv('/Users/tanishi/FYP/CDS/processed/NEG_4_isoforms.csv')

# ------------ Extract FASTA file with only the sequences of the longest isoforms 
def write_longest_isoforms_fasta(original_fasta_path, output_fasta_path, longest_isoform_ids):
    """
    Writes a new FASTA file containing only the sequences of the longest isoforms.

    Parameters:
    -----------
    original_fasta_path : str
        Path to the original FASTA file with all isoforms.

    output_fasta_path : str
        Path where the filtered FASTA file will be saved.

    longest_isoform_ids : set or list
        A collection of FBpp IDs (isoform IDs) that should be retained in the new FASTA.
    """
    with open(original_fasta_path, 'r') as infile, open(output_fasta_path, 'w') as outfile:
        write_seq = False
        for line in infile:
            if line.startswith('>'):
                current_id = line[1:].split()[0]
                write_seq = current_id in longest_isoform_ids
            if write_seq:
                outfile.write(line)

# Get the set of longest isoform IDs
longest_ids = set(NEG_4_isoforms['Longest_Isoform_ID'])

# Define input and output FASTA paths
input_fasta = "/Users/tanishi/FYP/CDS/NEG_CDS_4.fasta"
output_fasta = "/Users/tanishi/FYP/CDS/processed/NEG_isoform_4.fasta"

# Call the function
write_longest_isoforms_fasta(input_fasta, output_fasta, longest_ids)





