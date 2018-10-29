#! /usr/bin/env python3
import numpy as np
import pandas as pd

def read_mfasta(filename):
    """
    this function takes a mfasta file and returns sequences with all gaps "-"
    important note :only one line for each fasta sequence
    """
    seq_list = []
    with open (filename, 'r') as filin:
        for line in filin:
            if line[0] == ">":
                seq = filin.readline().rstrip()
                seq_list.append(seq)
    return seq_list

def count_gap(seq_list):
    """ this function creates a 1D array to count the gap and it returns a 1D array 
    containing the numbers of gaps of each column(position of a sequence) """
    gap_mt = np.zeros(len(seq_list[0]), int)
    cols_todel = []
    for seq in seq_list:
        for col_idx in range(len(seq_list[0])):
            if seq[col_idx] == "-":
                gap_mt[col_idx] += 1
    return gap_mt

def get_coltodel(mt):
    """get the column_index that the numbers of gaps is more than the half of length of aligmed seauences""" 
    cols_todel = []
    for pos in range(mt.shape[0]):
        count = mt[pos]
        if count >= 1/2*(len(seqList)):
            cols_todel.append(pos)
    return cols_todel

### main program

# create a list containing all fasta sequences with gaps
seqList = read_mfasta("query.v2_mfasta")

#create a dataframe with each row presenting a fasta sequences only with few gaps 
gaps = count_gap(seqList)
cols_todel = get_coltodel(gaps)
df = pd.DataFrame.from_records(seqList)
# remove the columns whose gaps are equal or more than half of numbers of input sequence
df.drop(cols_todel, axis =1, inplace=True)

# give the column names
col_names =[]
for i in range(1, (df.shape[1]+1)):
    name = "Pos"+ str(i)
    col_names.append(name)
df.columns = col_names
print("the obtained dataframe is :")
print(df)

# get the frequency of each aa per position and reindex the rows to have the same orders as aa_orders
aa_orders = "ACDEFGHIKLMNPQRSTVWY-"
count_table = df.apply(pd.value_counts).reindex(list(aa_orders)).fillna(0)
freq_table = count_table / count_table.sum(axis = 0)[0]
print("the obtained frequence table is :")
print(freq_table)

### next step 1: try to code a python script to get sequece_weight
### next step 2: figure out how to combine the seq_weight with freq_table
### next step 3: get background frequency of each amine acid