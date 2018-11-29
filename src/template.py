#! /usr/bin/env python3
import numpy as np
import pandas as pd
import sys

PATH_METAFOLD_LIST="/home/madeleine/Documents/2018---2019-partage/Data/METAFOLD_extended.list"
def get_name(filename):
    with open (filename, 'r') as f:
        line = f.readline()
        name = line[1:]
    return name

def read_map(filename):
    """
    this function takes template.map file and returns a list of sequences
    """
    chevrons=[]
    seq_list = []
    seq_list_name = []
    with open (filename, 'r') as f:
        lines = f.readlines()
        for num,line in enumerate(lines, 0):
            if line.startswith(">"):
                chevrons.append(num)
                seq_list_name.append(line.split(";")[-1].replace("\n",""))
    for c in range(0,len(chevrons) - 1 ):
            deb = chevrons[c]
            fin = chevrons[c+1]
            tmp=[]
            for x in range(deb + 2 ,fin ):
                tmp.append(lines[x].replace('\n', '').replace('*' ,''))
            seq_list.append(''.join(tmp))
    return(seq_list, seq_list_name)

def get_cols_todel(seq_list, index_real_query):

    """ this function performs 3 steps:
    1. creates a 1D array to count the gaps of each position of a given sequence
    2. creates a list to indicate the position SHOULD be keeped
    3. finally return a list of positions to delete when there are too many gaps in the position and this position is NOT keeped position
    (more than half of input sequences) """

    gap_mt = np.zeros(len(seq_list[0]), int) #create a 1D array for indicating all positions of sequence
    cols_todel = []
    for seq in seq_list:
        for col_idx in range(len(seq_list[0])):
            if seq[col_idx] == "-":
                gap_mt[col_idx] += 1

    #set a list to keep the positions that could not be deleted (the position in which one amine acid exists in sequence)
    seq_ref = seq_list[index_real_query]
    print(seq_ref)
    cols_keep = []
    for seq_pos in range(len(seq_ref)):
        if seq_ref[seq_pos] != "-":
            cols_keep.append(seq_pos)

    #get the column_index that the numbers of gaps are more than half of length of input sequences and NOT in the positions to keep
    for pos in range(gap_mt.shape[0]):
        count = gap_mt[pos]
        if count >= 1/2*(len(seq_list)) and pos not in cols_keep:
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


def create_pssm(seq_list, wgt_list):

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
            array[pos, aa_idx] += round(1*seq_wgt, 5)

    # add the pseudo-count (background frequency of aa = 1/20) in each column
    bg = np.full((len(seq_list[0]), 20), 0.05000)
    gap = np.zeros((len(seq_list[0]), 1))
    array_bg = np.append(bg, gap, axis = 1)
    new_array = (array + array_bg)/2
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


def find_who_is_the_real_seq(template_name):
    #### A SUPPRIMER QUAND VOUS AUREZ A FAIRE TOURNER LE SCRIPT (j'ai changé le
    ### nom du template parce que je suis dans un dossier test) ####
    template_name=template_name.split('/')[-1].split(".")[0]
    ############################################################################
    with open(PATH_METAFOLD_LIST, "r") as fillin:
        for elem in fillin:
            if elem.split()[0] == template_name :
                return elem.split()[1].split(".")[0]

if __name__ == "__main__":

    aa_orders = "ARNDCQEGHILKMFPSTWYV-"

    # create a list containing all fasta sequences with gaps
    seqList, seqListName = read_map(sys.argv[1])
    template_name = sys.argv[1]

    #find where is our query in .map
    real_query_name_in_map = find_who_is_the_real_seq(template_name)
    index_query_name_in_map = seqListName.index(real_query_name_in_map)

    #create a dataframe with each row presenting a fasta sequences only with few gaps
    pos_todel = get_cols_todel(seqList, index_query_name_in_map)
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

    # get the weights of each conserved sequence
    seq_wgts = get_weights(conserved_seq, freq_list)

    # get the position where there is a gap in the onserved sequence
    pos_out_pssm = []
    cons_seq = conserved_seq[index_query_name_in_map]
    for pos_idx in range(len(cons_seq)):
        if cons_seq[pos_idx] == "-":
            pos_out_pssm.append(pos_idx)
    cons_seq= cons_seq.replace("-","")

    # create the pssm
    M_pssm = create_pssm(conserved_seq, seq_wgts)
    pssm = pd.DataFrame(M_pssm)
    pssm.columns = list(aa_orders)

    # check if there is a gap in query sequence, we remove this position in pssm
    if len(pos_out_pssm) != 0:
        new_pssm = pssm.drop(pssm.index[pos_out_pssm])
        new_pssm = new_pssm.reset_index(drop=True)
    else:
        new_pssm=pssm
    ss2 = get_ss2_df(sys.argv[2])
    frames = [new_pssm, ss2]
    pssm_ss2 = pd.concat(frames, axis=1, ignore_index=True)

    with open(sys.argv[3], "w") as output_f:
        print(template_name)
        output_f.write(template_name.split(".")[0]+"\n")
        output_f.write(cons_seq +"\n")
        pssm_ss2.to_string(output_f, header=False, index=False)
