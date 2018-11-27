#! /usr/bin/env python3
from scoring_matrix import ScoringMatrix
import numpy as np
import terminal_output
import glob
import sys
import os
from align_pssm import *
import random as rd
import numpy as np

gap_score = -1
terminal_gap_score = 0
gap_open = - 0.5
blossum = "/home/sdv/m2bi/stieo/2018---2019-Equipe-2/data/data_test/blossum62.txt"

def blossum_scores_val( blossum ):
    with open (blossum, 'r') as f:
        scores = f.read()
        scores = list(map(int, scores.split()))
    blossum_mean = np.mean(scores)
    blossum_std = np.std(scores)
    A = (10 - blossum_mean)/blossum_std
    return(A)

def gen_random_scores (pssm1,pssm2):
    random_scores = []
    for i in range (0,10000):
        r1 = rd.randint(0,len(pssm1[2]) -1 )
        r2 = rd.randint(0,len(pssm2[2]) -1 )
        random_scores.append(score_between_2_pssm(pssm1[2][r1],pssm2[2][r2]))
    return(random_scores)


def optimal_gap ( pssm1, pssm2):
    random_scores = gen_random_scores (pssm1,pssm2)
    A = blossum_scores_val( blossum )
    gap_open = ( A + np.mean(random_scores) ) * np.std(random_scores)
    gap_extension = (gap_open * 10)/100
    return([gap_open, gap_extension])


if __name__ == "__main__":
    cwd = os.getcwd()
    path="/home/sdv/m2bi/stieo/2018---2019-Equipe-2/data/pssm_templates/"
    all_open_gap=[]
    all_ext_gap =[]
    for i in range(0,100):
        path_pssm1 = path + rd.choice(os.listdir("/home/sdv/m2bi/stieo/2018---2019-Equipe-2/data/pssm_templates/"))
        path_pssm2 = path + rd.choice(os.listdir("/home/sdv/m2bi/stieo/2018---2019-Equipe-2/data/pssm_templates/"))
        pssm1= read_pssm(path_pssm1)
        pssm2=read_pssm(path_pssm2)
        print(optimal_gap (pssm1, pssm2))
        all_open_gap.append(optimal_gap (pssm1, pssm2)[0])
        all_ext_gap.append(optimal_gap (pssm1, pssm2)[1])
    print(np.mean(all_open_gap))
    print(np.std(all_ext_gap))
