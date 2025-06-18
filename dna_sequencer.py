# import Bio
from Bio.Seq import Seq
# from Bio.SeqUtils import GC
# The function GC was part of the Bio.SeqUtils module in older versions of Biopython. 
# It has been deprecated and removed in favor of a more intuitive, object-oriented approach.
# Instead of calling a separate function, you now call the .gc_content() method directly on 
# a Seq object.
from Bio.SeqUtils import gc_fraction
from collections import Counter
from Bio import SeqIO

def validate_dna_sequence(sequence):
    # check if sequence is DNA (only uppercase A, C, G, T)
    valid = set('ACGT')
    if set(str(sequence).upper()) != valid:
    #condition below is same as condition above
    # if not set(str(sequence).upper()).issubset(valid):
        raise ValueError("Invalid DNA sequence. Only A, G, C, T allowed.")
    return sequence.upper()

def read_dna_sequence(file_path):
    with open(file_path, 'r') as file:
        sequence = file.read().strip()
    return sequence

def main(file_path):
    sequence = read_dna_sequence(file_path)
    validate_dna_sequence(sequence)
    print(f"Sequence: {sequence}")

if __name__ == "__main__":

    file_path = "sequence.txt"

    main(file_path)
