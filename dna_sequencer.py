# import Bio
from Bio.Seq import Seq
# from Bio.SeqUtils import GC
# The function GC was part of the Bio.SeqUtils module in older versions of Biopython. 
# It has been deprecated and removed in favor of a more intuitive, object-oriented approach.
# Instead of calling a separate function, you now call the .gc_content() method directly on 
# a Seq object.
from Bio import SeqIO

def read_dna_sequence(file_path):
    with open(file_path, 'r') as file:
        sequence = file.read().strip()
    return sequence

def main(file_path):
    sequence = read_dna_sequence(file_path)
    print(f"Sequence: {sequence}")

if __name__ == "__main__":

    file_path = "sequence.txt"

    main(file_path)
