from Bio.Seq import Seq

def reverse_complement(sequence):
    biopython_seq = Seq(sequence)
    return str(biopython_seq.reverse_complement())
