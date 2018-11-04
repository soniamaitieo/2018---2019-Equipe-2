#! /usr/bin/env python3
import numpy as np
import os


def score_between_2_pssm(v1,v2,mblos_path):
    """
    this function goal is to calculate score bewteen 2 vectors

    """
    score=0
    m_blos=np.loadtxt(mblos_path)

    for i in range(0,len(v1)) :
        for j in range(0,len(v2)):
            score=score+v1[i]*v2[i]*m_blos[i][i]
    return score

def dynamic_alignment (pssm1,pssm2, gap_score,mblos_path) :
    n=pssm1.shape[0]
    m=pssm2.shape[0]
    print(n,m)
    align_matrix=np.zeros((n,m), dtype=float)
    align_matrix[0,0]=0
    for i in range(1,n) : align_matrix[i,0]=0
        #align_matrix[i,0]=align_matrix[i-1,0]+gap_score
    for j in range(1,m) : align_matrix[0,j]=0
        #align_matrix[0,j]=align_matrix[0,j-1]+gap_score

    for i in range(1,n) :
        for j in range(1,m) :
            match=align_matrix[i-1,j-1]+score_between_2_pssm(pssm1[i],pssm2[j],mblos_path)
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



def graph() :
    from statistics import mean
    from statistics import var
    all_DA = choose_gap_op(PSSM1,PSSM2)

    gap1 = [j for i in all_DA[0].tolist() for j in i]
    gap2 = [j for i in all_DA[1].tolist() for j in i]
    gap3 = [j for i in all_DA[2].tolist() for j in i]
    gap4 = [j for i in all_DA[3].tolist() for j in i]
    gap5 = [j for i in all_DA[4].tolist() for j in i]

    all_scores_gaps=[gap1,gap2,gap3,gap4,gap5]
    all_scores_gaps_list=[j for i in truc for j in i]

    # matplotlib histogram : distribution des scores
    plt.hist(all_scores_gaps_list , color = 'blue', edgecolor = 'black',
         bins = int(180/5))
    mean(trucc)
    np.var(trucc)



if __name__ == "__main__":
    cwd = os.getcwd()
    name_fasta1,seq1,PSSM1=read_pssm(cwd+"/data/query_small.aamtx")
    name_fasta2,seq2,PSSM2=read_pssm(cwd+"/data/template_small.aamtx")
    #score_between_2_pssm(PSSM1,PSSM2)
    mblos_path=cwd+"/data/BLOSUM62.txt"
    DA = dynamic_alignment(PSSM1,PSSM2, -1,mblos_path)
    print(DA)
    #get_alignments(DA,seq1,seq2)
    #caca = dynamic_alignment(PSSM1,PSSM2, -2)



######################a revoir ############################################

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


    #####################################################################
