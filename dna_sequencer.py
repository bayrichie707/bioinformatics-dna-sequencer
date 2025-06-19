# import Bio
import argparse
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

def analyze_dna_sequence(file_path):
    # Analysis function using biopython
    try:
        # Use SeqIO to support formats such as FASTA
        # Use parse to check file for multiple lines
        for record in SeqIO.parse(file_path, "fasta"):
            sequence = validate_dna_sequence(record.seq)

        stats = {
            'Sequence ID': record.id,
            'Sequence Length': len(sequence),
            'GC Ratio': gc_fraction(sequence) * 100
        }

        # Call find_motifs of specific length. Return dictionary of motifs with more than one occurence.
        motifs = find_motifs(sequence, 3)
        stats['Motifs (length = 3)'] = motifs

        return stats
    
    except FileNotFoundError:
        print(f"Error: File '{file_path} not found.")

def find_motifs(sequence, motif_length):
    # Count occurrenes of motifs
    motifs = {}
    for i in range(len(sequence) - motif_length + 1):
        # Build the dictionary and store all substrings of length 4 from the sequence 
        motif = str(sequence[i : i + motif_length])
        # Count how many times motif (ATG, TGC) appears. If 0 then increment by 1
        motifs[motif] = motifs.get(motif, 0) + 1
    
    # Only return motifs that are greater than one
    return {motif: count for motif, count in motifs.items() if count > 1}




def read_dna_sequence(file_path):
    with open(file_path, 'r') as file:
        sequence = file.read().strip()
    return sequence

def print_results(stats):
    print("DNA SEQUENCE ANALYSIS RESULTS")
    print(f"Sequence ID: {stats['Sequence ID']}")
    print(f"Sequence Length: {stats['Sequence Length']}")
    print(f"GC Ratio: {stats['GC Ratio']:.2f}%")
    print("***Motifs (length = 3)***")
    motifs = stats['Motifs (length = 3)']
    if motifs:
        for motif, count in motifs.items():
            print(f"  {motif}: {count} times")
    else:
        print("  None found.")

    # print(motifs)

# def main(file_path):

    # sequence = read_dna_sequence(file_path)
    # validate_dna_sequence(sequence)
    # analyze_dna_sequence(file_path)
    # print(f"Sequence: {sequence}")
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='DNA Sequencer with Biopython')
    parser.add_argument('file', help='Path to the DNA Seq file (FASTA format)')
    args = parser.parse_args()
    
    stats = analyze_dna_sequence(args.file)
    # file_path = "sequence.txt"
    # file_path = "chromosome 11.fasta"
    print_results(stats)
    #main(file_path)
