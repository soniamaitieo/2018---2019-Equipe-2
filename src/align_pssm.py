#! /usr/bin/env python3
import numpy as np

# mt = np.zeros((len(seq_list[0]), 20), int)

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

def score_between_2_pssm(pssm1,pssm2):
    """
    this function goal is to look at all collum in 2 pssm
    and look at the correlation score of the 2 pssm

    """
    col_pssm1=[]
    col_pssm2=[]
    score=0
    for i in range(0,len(pssm1[1])) :
        col_pssm1.append(pssm1[1][i])
        col_pssm2.append(pssm2[2][i])
        #analyse de score ici



if __name__ == "__main__":
    seq_list = read_mfasta("query.v2_mfasta")
    mt1 = np.zeros((len(seq_list[0]), 20), int)
    aa_column = "ACDEFGHIKLMNPQRSTVWY"
    for seq in seq_list:
        for pos in range(len(seq)):
            aa_idx = aa_column.find(seq[pos])
            mt1[pos, aa_idx] += 1
    mt2 = np.zeros((len(seq_list[0]), 20), int)
    score_between_2_pssm(mt1,mt2)
