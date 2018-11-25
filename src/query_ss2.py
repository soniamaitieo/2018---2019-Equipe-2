#! /usr/bin/env python3
import numpy as np
import pandas as pd
import sys

def get_query_name(filename):
    with open (filename, 'r') as f:
        line = f.readline()
        query_nqme = line[1:]
    return query_nqme

def read_mfasta(filename):
    """
    this function takes a mfasta file and returns sequences with all gaps "-"
    important note :only one line for each fasta sequence (\n has been removed)
    """
    seq_list = []
    with open (filename, 'r') as filin:
        for line in filin:
            if line[0] == ">":
                seq = filin.readline().rstrip()
                seq_list.append(seq)
    return seq_list

def get_cols_todel(seq_list):
    """ this function firstly creates a 1D array to count the gaps of each position of a sequence.
    it returns subsequently a list of positions to delete when there are too many gaps in one position
    (more than half of input sequences) """
    gap_mt = np.zeros(len(seq_list[0]), int)
    cols_todel = []
    for seq in seq_list:
        for col_idx in range(len(seq_list[0])):
            if seq[col_idx] == "-":
                gap_mt[col_idx] += 1

    #get the column_index that the numbers of gaps are more than half of length of input sequences
    for pos in range(gap_mt.shape[0]):
        count = gap_mt[pos]
        if count >= 1/2*(len(seq_list)):
            cols_todel.append(pos)
    return cols_todel

def get_conserved_seq(table):
    """ create a dictionary to get the conserved sequences """
    seq_dict = {}
    for idx, row in table.iterrows():
        seq = ""
        for pos in table.columns:
            seq += row[pos]
        seq_dict[idx] = seq
    return seq_dict

def get_freq(table):
    """ to get the frequency of each columns (positions) and put into a list"""
    pos_freq_lst = []
    for pos in table.columns:
        sequence = "".join(table[pos].values)
        freq = {aa: sequence.count(aa) for aa in sequence}
        pos_freq_lst.append(freq)
    return pos_freq_lst

def get_conserved_seq(table):
    """ create a list to get the conserved sequences """
    seq_list = []
    for idx, row in table.iterrows():
        seq = ""
        for pos in table.columns:
            seq += row[pos]
        seq_list.append(seq)
    return seq_list

def get_weights(seq_list, freq_list):
    """ takes 2 lists (list of sequences and list of aa frequency of each position) as input
        returns a list of weights of input sequences"""
    wgt_list = []
    for sequ in seq_list:
        total = 0
        count = 0
        for id in range(len(sequ)):
            for key, value in freq_list[id].items():
                if sequ[id] == key:
                    prop = 1/(value*len(freq_list[id]))
                    total += prop
                    count += 1
        weight = round(total/count, 5)
        wgt_list.append(weight)
    return wgt_list


def create_pssm(seq_list, wgt_list, bg_freq):

    """ this function creates firstly a matrix of nx21 (n is the length of input seq)
    each input sequence is compared to aa_orders = "ARNDCQEGHILKMFPSTWYV-" to calculate
    the frequences of aa in each position
    """
    array = np.zeros((len(seq_list[0]), 21))
    for idx in range(len(seq_list)):
        seq = seq_list[idx]
        seq_wgt = wgt_list[idx]
        for pos in range(len(seq)):
            aa_idx = aa_orders.find(seq[pos])
            array[pos, aa_idx] += 1*seq_wgt

    # add the pseudo-count (background frequency of aa) in each column
    for order in range(len(aa_orders[:-1])):
        for aa, freq in bg_freq.items():
            if aa_orders[order] == aa:
                array[:, order] += round(bg_freq[aa], 5)
    new_array = array/2
    return new_array

def get_ss2_df(ss2_file):
    # reads output file returned by psipred and returns output.ss2 by removing non-useful information
    with open (ss2_file, 'r') as filin, open("output.ss2", "w") as filout:
        for lines in filin.readlines():
            if len(lines) > 1 and lines.strip()[0] != "#":
                filout.write(lines.strip()+"\n")

    # create a dataframe by taking output.ss2 as input
    ss2_df = pd.read_table("output.ss2", delim_whitespace=True, header=None)
    ss2_df.drop([0,1,2], axis =1, inplace=True) #remove the first 3 columns
    ss2_df.columns = ["C", "H", "E=strand"] # give the names of column corresponding to secondary structure
    return ss2_df


if __name__ == "__main__":

    aa_orders = "ARNDCQEGHILKMFPSTWYV-"

    bg_freq = {"A":0.0789, "R":0.054, "N": 0.0413, "D": 0.0535, "C": 0.0150, "Q":0.0395, "E": 0.0667,
               "G":0.0696, "H":0.0229, "I": 0.059, "L": 0.0965, "K": 0.0592, "M":0.0238, "F": 0.0396,
               "P":0.0483, "S":0.0682, "T": 0.0541, "W": 0.0113, "Y":0.0303, "V":0.0673}

    # create a list containing all fasta sequences with gaps

    seqList = read_mfasta(sys.argv[1])

    #create a dataframe with each row presenting a fasta sequences only with few gaps
    pos_todel = get_cols_todel(seqList)
    df = pd.DataFrame.from_records(seqList)
    # remove the columns whose gaps are equal or more than half of numbers of input sequence
    df.drop(pos_todel, axis =1, inplace=True)

    # give the column names
    col_names =[]
    for i in range(1, (df.shape[1]+1)):
        name = "Pos"+ str(i)
        col_names.append(name)
    df.columns = col_names

    # get the frequency of each columns (positions)
    freq_list = get_freq(df)

    # get the conserved sequences of mfasta in which most gaps are removed
    conserved_seq = get_conserved_seq(df)

    # get the weights of each conserved sequences
    seq_wgts = get_weights(conserved_seq, freq_list)

    # create the pssm
    M_pssm = create_pssm(conserved_seq, seq_wgts, bg_freq)
    pssm = pd.DataFrame(M_pssm)
    pssm.columns = list(aa_orders)
    ss2 = get_ss2_df(sys.argv[2])
    frames = [pssm, ss2]
    pssm_ss2 = pd.concat(frames, axis=1)

    with open(sys.argv[3], "w") as output_f:
        output_f.write(get_query_name(sys.argv[1]))
        output_f.write(conserved_seq[0]+"\n")
        pssm_ss2.to_string(output_f, header=False, index=False)
