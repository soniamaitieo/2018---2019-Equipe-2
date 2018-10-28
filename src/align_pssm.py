#! /usr/bin/env python3
import numpy as np

def score_between_2_pssm(v1,v2):
    """
    this function goal is to calculate score bewteen 2 vectors

    """
    score=0
    m_blos=np.loadtxt("/home/soniamai/Bureau/2018---2019-Equipe-2/data/BLOSUM62.txt")

    for i in range(0,len(v1)) :
        for j in range(0,len(v2)):
            score=score+v1[i]*v2[i]*m_blos[i][i]
    return score

def dynamic_alignment (pssm1,pssm2, gap_score) :
    n=pssm1.shape[0]
    m=pssm2.shape[0]
    print(n,m)
    align_matrix=np.zeros((n,m), dtype=float)
    align_matrix[0,0]=0
    for i in range(1,n) : align_matrix[i,0]=align_matrix[i-1,0]+gap_score
    for j in range(1,m) : align_matrix[0,j]=align_matrix[0,j-1]+gap_score

    for i in range(1,n) :
        for j in range(1,m) :
            match=align_matrix[i-1,j-1]+score_between_2_pssm(pssm1[i],pssm2[j])
            align_matrix[i,j]=max(match,align_matrix[i-1,j]+gap_score,align_matrix[i,j-1]+gap_score)
    return(align_matrix)
    
def read_pssm(path) :
    ##1ere ligne : nom du fasta
    ##2eme ligne : sequence
    ##a partir de la 3eme : tableau
    with open(path,"r") as fillin :
        name_fasta = fillin.readline()
        seq=fillin.readline()
        lines = (line for line in fillin)
        FH = np.loadtxt(lines)
    return(name_fasta,seq,FH)
    
    
def choose_gap_op(pssm1,pssm2):
    L=[]
    gaps=[-1,-2,-3,-4,-5]
    for gap in gaps :
        L.append(dynamic_alignment (pssm1,pssm2, gap))
    return(L)

all_DA = choose_gap_op(pssm1,pssm2)


# Import the libraries
import matplotlib.pyplot as plt

# matplotlib histogram
plt.hist(list(all_DA[1][3]), color = 'blue', edgecolor = 'black',
         bins = int(180/5))



 
def score_shuffle_pssm(PSSM1,PSSM2) :
    pssm1shuffle  = np.array(PSSM1)
    pssm2shuffle =  np.array(PSSM2)
    np.random.shuffle(pssm1shuffle)
    np.random.shuffle(np.array(PSSM2))
    random_da = dynamic_alignment(pssm1shuffle, pssm2shuffle)
    return(random_da[PSSM1.shape[0] -1,PSSM2.shape[0]-1])

def n_random_scores(n) :
    for i in range (1,n):
        print(score_shuffle_pssm(PSSM1,PSSM2))
    
    

    
if __name__ == "__main__":
    name_fasta1,seq1,PSSM1=read_pssm("/home/soniamai/Bureau/2018---2019-Equipe-2/data/query_small.aamtx")
    name_fasta2,seq2,PSSM2=read_pssm("/home/soniamai/Bureau/2018---2019-Equipe-2/data/template_small.aamtx")
    #score_between_2_pssm(PSSM1,PSSM2)
    DA = dynamic_alignment(PSSM1,PSSM2, -2)
