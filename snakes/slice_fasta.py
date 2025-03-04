#!/usr/bin/env python

import sys

try:
    infilename = sys.argv[1]
    outfilepref = sys.argv[2]
except:
    print("""
Repartir secuencias de un multifasta en 10 rebanadas.
Ejemplo:
slice_fasta.py <infile.fastq.gz> <outfile_prefix>
    """)
    sys.exit()

def divide_multifasta(input_file, output_prefix, num_slices=10):
    """Divide a multifasta file into a specified number of slices."""
    
    # Read the entire multifasta file
    with open(input_file, 'r') as fasta_file:
        fasta_contents = fasta_file.read()

    # Split the contents into individual sequences based on '>' symbol
    sequences = fasta_contents.strip().split('>')
    
    # Remove any empty entries (this can happen if the file starts with '>')
    sequences = [seq for seq in sequences if seq]

    # Calculate the number of sequences per slice
    total_sequences = len(sequences)
    sequences_per_slice = (total_sequences + num_slices - 1) // num_slices  # Ceiling division

    # Create the slices
    for i in range(num_slices):
        # Calculate start and end index for this slice
        start_index = i * sequences_per_slice
        end_index = min(start_index + sequences_per_slice, total_sequences)

        # If start_index is beyond the total, we have no more sequences to write
        if start_index >= total_sequences:
            break

        # Prepare the output filename
        output_file = f"{output_prefix}_slice_{i + 1}.fasta"

        # Write the slice to the output file
        with open(output_file, 'w') as slice_file:
            for sequence in sequences[start_index:end_index]:
                slice_file.write(f'>{sequence.strip()}\n')  # Write sequence back with '>' symbol

        print(f'Created {output_file} with sequences {start_index + 1} to {end_index}')

# Example usage
if __name__ == "__main__":
    divide_multifasta(infilename, outfilepref)