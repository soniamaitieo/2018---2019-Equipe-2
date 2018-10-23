#! /usr/bin/env python3
import numpy as np

def read_mfasta(filename):
    """
    this function takes a mfasta file and returns sequences without gaps "-"
    notes :
    1. only one line for each fasta sequence
    2. for the moment just output the sequences with the same length of query >>> needed to improve the codes
    """
    seq_list = []
    with open (filename, 'r') as filin:
        for line in filin:
            if line[0] == ">" and line.find("73 aa") != -1:
                seq = filin.readline().rstrip()
                clean_seq = seq.replace("-", "")
                seq_list.append(clean_seq)
    return seq_list

def mt_of_aa_freq(seq_list):

    """ this function creates firstly a matrix of nx20 (n is the length of input seq)
    each input sequence is compared to aa_column = "ACDEFGHIKLMNPQRSTVWY" to calculate
    the frequences of aa in each position
    """
    mt = np.zeros((len(seq_list[0]), 20), int)
    for seq in seq_list:
        for pos in range(len(seq)):
            aa_idx = aa_column.find(seq[pos])
            mt[pos, aa_idx] += 1
    return mt


if __name__ == "__main__":
    
    aa_column = "ACDEFGHIKLMNPQRSTVWY"

    mylist = read_mfasta("query.v2_mfasta")
    print(len(mylist))
    print(len(mylist[1]))
    
    mt_freq = mt_of_aa_freq(mylist)
    print(mt_freq.shape)
    print(mt_freq)