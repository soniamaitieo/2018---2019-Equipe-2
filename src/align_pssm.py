#! /usr/bin/env python3

# This file is part of sequence-aligner.
# Copyright (C) 2014 Christopher Kyle Horton <chorton@ltu.edu>
# Modified by Madeleine DE SOUSA


import os
import sys
import glob
import numpy as np
from scoring_matrix import ScoringMatrix


GAP_SCORE = -1
GAP_OPEN = - 0.5
TERMINAL_GAP_SCORE = 0

CWD = os.getcwd()
HOMSTRAD_PATH = "/home/madeleine/Documents/2018---2019-partage/Data/HOMSTRAD/"

def score_between_2_pssm(v1, v2):
    """
    this function goal is to calculate score bewteen 2 vectors

    """
    score = 0
    for i in range(len(v1)):
        score = score+v1[i]*v2[i]
    return score

def score_Pearson(v1, v2, pssm1, pssm2):
    """
    this function goal is to calculate score between 2 vector with Pearsons's correlation
    """
    score_up = 0
    score_down_1 = 0
    score_down_2 = 0
    moy_pssm1 = np.mean(pssm1, axis=0)
    moy_pssm2 = np.mean(pssm2, axis=0)
    for i in range(len(v1)):
        score_up += (v1[i]-moy_pssm1[i]) * (v2[i]-moy_pssm2[i])
        score_down_1 += (v1[i]-moy_pssm1[i])**2
        score_down_2 += (v2[i]-moy_pssm2[i])**2
    return score_up/np.sqrt(score_down_1*score_down_2)

def initialize_edges(sm):
    '''Sets up the top and left edges of the provided ScoringMatrix. This is
    the first step in the dynamic programming algorithm.
    Performs a semi-global alignment by default unless alignment_is_global is
    specified to be True.'''
    sm.set_score(0, 0, 0)
    gap = TERMINAL_GAP_SCORE
    for i in range(1, sm.get_rows()):
        sm.set_score(i, 0, sm.get_score(i - 1, 0) + gap)
        sm.add_up_backlink(i, 0)
    for i in range(1, sm.get_columns()):
        sm.set_score(0, i, sm.get_score(0, i - 1) + gap)
        sm.add_left_backlink(0, i)

def fill_matrix(sm, pssm1, pssm2):
    '''Uses dynamic programming to fill out a provided ScoringMatrix after
    the edges have been initialized.
    Performs a semi-global alignment by default unless alignment_is_global is
    specified to be True.
    Based on pseudocode provided on page 54 of our textbook.'''
    initialize_edges(sm)
    score_diagonal = 0

    for i in range(1, sm.get_rows() +1):
        for j in range(1, sm.get_columns() +1):
            # Calculate scores
#            match=score_between_2_pssm(pssm1[i-1],pssm2[j-1])
            match = score_Pearson(pssm1[i - 1], pssm2[j - 1], pssm1, pssm2)
#            print("score match dot_product: " + str(match))
#            print("score match pearsons : " + str(match2))
            score_diagonal = sm.get_score(i - 1, j - 1) + match
            if sm.get_backlinks(i, j - 1)["left"]:
                score_left = sm.get_score(i, j - 1) + GAP_OPEN
            else:
                score_left = sm.get_score(i, j - 1) + GAP_SCORE

            if sm.get_backlinks(i - 1, j)["up"]:
                score_up = sm.get_score(i - 1, j) + GAP_OPEN
            else:
                score_up = sm.get_score(i - 1, j) + GAP_SCORE


            max_score = max(score_diagonal, score_left, score_up)
            sm.set_score(i, j, max_score)
            # Establish backlink(s)
            if max_score == score_diagonal:
                #print(match)
                sm.add_diagonal_backlink(i, j)
            if max_score == score_left:
                sm.add_left_backlink(i, j)
            if max_score == score_up:
                sm.add_up_backlink(i, j)

def get_alignments(sm, pssm1, pssm2):
    '''Returns a list of the alignments generated from the scoring matrix.
    Performs a semi-global alignment by default unless alignment_is_global is
    specified to be True.'''
    fill_matrix(sm, pssm1, pssm2)
    # Find the optimal alignment paths starting from the lower right corner.
    seq = (sm.get_top_sequence(), sm.get_left_sequence())
    last_col = len(seq[0])
    last_row = len(seq[1])

    todo_list = [[last_row, last_col, "", ""]] # Entry (row,col,string0,string1)
    done_list = [] # Entry (str0,str1)
    backlink_used = [[False for x in range(len(seq[0]) + 1)] for x in range(len(seq[1]) + 1)]
    while todo_list:
        row, col, str0, str1 = todo_list.pop()
        backlinks = sm.get_backlinks(row, col)
        if True in backlinks.values(): # If some back-link exists.
            backlink_used[row][col] = True # Mark linked cells as used as we go.
            if backlinks["diagonal"]:
                todo_list.append([row - 1, col - 1, seq[0][col - 1] + str0, seq[1][row - 1] + str1])
            if backlinks["up"]:
                todo_list.append([row - 1, col, '-' + str0, seq[1][row - 1] + str1])
            if backlinks["left"]:
                todo_list.append([row, col - 1, seq[0][col - 1] + str0, '-' + str1])
        else:
            done_list.append([str0, str1])
    # Clean up unused backlinks
    for row in range(len(seq[1]) + 1):
        for col in range(len(seq[0]) + 1):
            if not backlink_used[row][col]:
                sm.remove_backlinks(row, col)

    return (done_list, sm.get_score(row, col))

def read_pssm(path):
    """
    Function which read PSSM from a file
    """
    ##1ere ligne : nom du fasta
    ##2eme ligne : sequence
    ##a partir de la 3eme : tableau
    with open(path, "r") as fillin:
        name_fasta = fillin.readline()
        seq = fillin.readline().replace(" ", "").replace("\n", "")
        lines = (line for line in fillin)
        matrix = np.loadtxt(lines)
    #print(seq)
    list_a_suppr = []
    #for i in range(0,len(seq)):
    for i, sequence in enumerate(seq):
        if sequence == "-":
            list_a_suppr.append(i)
    #print(FH.shape)
    #print(list_a_suppr)
    for elem in reversed(list_a_suppr):
        matrix = np.delete(matrix, elem, axis=0)
    #print(FH.shape)
    return(name_fasta, seq.replace("-", ""), matrix)

def psipred_for_foldrec(pssm):
    list_score_psipred = []
    list_struc_psipred = []
    for line in pssm:
        coil = line[-3]
        helix = line[-2]
        sheet = line[-1]

        if coil >= helix and coil >= sheet:
            list_struc_psipred.append("C")
        elif helix >= coil and helix >= sheet:
            list_struc_psipred.append("H")
        else:
            list_struc_psipred.append("E")
        list_score_psipred.append(max(coil, helix, sheet))
    return((list_struc_psipred, list_score_psipred))

def make_output(name_query, list_score, dico_all, seq1, dico_psipred):
    """
    Prepare output_file with all alignments prepared before and their score
    """


    with open(name_query+".foldrec1", "w") as fillout1:
        with open(name_query+".foldrec2", "w") as fillout2:
            fillout1.write("*** HITS RANKED *** \n\n")
            fillout1.write("SEQUENCE QUERY FILE : " + name_query + ", "+ str(len(seq1))+"\n\n")
            fillout1.write("#| Score | Ungaped_score | Pvalue_Q | Pscore | PQTscore | P-Value_T  \
            | Q. Length | T. Length | Q. begin-end | T. begin-end | HITS\n")
            fillout1.write("---------------------------------------------------------\n")

            fillout2.write("*** ALIGNMENTS DETAILS ***\n")

            number_template = 1
            struct_query, score_psipred = dico_psipred[name_query]

            for score in list_score:
                seq = dico_all[score][0][0]
                seq = seq.pop()
                name_template = dico_all[score][0][1]

                psipred_query = ""
                score_psipred_query = ""
                struct_query_tmp = struct_query.copy()
                score_psipred_tmp = score_psipred.copy()
                for elem in seq[1]:
                    if elem == "-":
                        psipred_query += "-"
                        score_psipred_query += "-"
                    else :
                        psipred_query += struct_query_tmp.pop(0)
                        score_psipred_query += str(score_psipred_tmp.pop(0))

                struct_temp, score_temp = dico_psipred[name_template]
                psipred_template = ""
                score_psipred_template = ""
                for elem in seq[0]:
                    if elem == "-":
                        psipred_template += "-"
                        score_psipred_template += "-"
                    else :
                        psipred_template += struct_temp.pop(0)
                        score_psipred_template += str(score_temp.pop(0))


                #with open(HOMSTRAD_PATH+name_template+"/"+name_template+".scop_id") as f :
                #    scop_id=f.readline().replace("\n","")


                fillout1.write("{:5}{:9}{:>8}{:>9}{:>8}{:>8}{:>6}{:>8}{:>6}{:>6}{:>9} \
                {:>9}{:^18}\n".format(str(number_template), round(score, 4), "X", "X", \
                 "X", "X", "X", "X", "X", "X", "X", "X", name_template))

                fillout2.write("No " + str(number_template) + "\n")

                fillout2.write("Alignment : " + name_query + ", " + str(len(seq1)) + \
                 " aa. vs "+ name_template + "\n")

                fillout2.write("Score : {:^10} | Normalized score : {:^10} | \
                 Query coverage : {:^10} | Identity : {:^10} | Gaps : {:^10}  \
                 | SS Score : {:10} |  Alignment length : {:10} | Corr score : {:10} \n\n" \
                 .format(score, "X", "X", "X", "X", "X", "X", "X"))

                fillout2.write("{:10} {:3} {:} {:>8}\n".format("Query", "1", seq[1], \
                 str(len(seq1.replace("_", "").replace("-", "")))))

                fillout2.write("{:10} {:3} {:} {:>8}\n".format("Query", "1", psipred_query, \
                 str(len(seq1.replace("_", "").replace("-", "")))))

                fillout2.write("{:10} {:3} {:} {:>8}\n".format("Query", "1", score_psipred_query,  \
                 str(len(seq1.replace("_", "").replace("-", "")))))

                fillout2.write("\n")

                fillout2.write("{:10} {:3} {:} {:>8}\n".format("Template", "1", seq[0], \
                 str(len(seq[0].replace("_", "")))))

                fillout2.write("{:10} {:3} {:} {:>8}\n".format("Template", "1", psipred_template, \
                 str(len(seq[0].replace("_", "")))))

                fillout2.write("{:10} {:3} {:} {:>8}\n".format("Template", "1", score_psipred_template, \
                 str(len(seq[0].replace("_", "")))))

                fillout2.write("\n")

                fillout2.write("{:10}\n".format("X " * len(seq[0])))
                fillout2.write("\n")

                del dico_all[score][0]
                number_template = number_template + 1
            fillout1.write("\n\n\n")

def run_alignment_and_make_output(fasta_file):
    """
    Main script which prepare a scoring matrix for all templates, run semi_global
    alignement and make output_file
    """
    name_fasta1, seq1, pssm1 = read_pssm(CWD + "/" + fasta_file)
    name_query = name_fasta1.replace(">", "").replace("\n", "").replace("'", "").replace(" ", "")
    dico_all = {}
    list_score = []
    dico_psipred = {}
    iterator = 1
    #CWD+"/data/data_test/rep_pssm_test/test/*aatmx"
    #CWD+"/data/pssm_templates/*aatmx"
    dico_psipred[name_query] = psipred_for_foldrec(pssm1)
    for elem in glob.glob(CWD + "/data/data_test/rep_pssm_test/test/*aatmx"):
        print("template number :" + str(iterator))
        iterator += 1
        name_fasta2, seq2, pssm2 = read_pssm(elem)
        name_fasta2 = name_fasta2.replace(">", "").replace("\n", "").replace("'", "") \
         .replace(" ", "").split(";")[-1]
        dico_psipred[name_fasta2]=psipred_for_foldrec(pssm2)


        matrix_scoring = ScoringMatrix(seq1, seq2)
        align, score = get_alignments(matrix_scoring, pssm1, pssm2)
        if score not in dico_all:
            dico_all[score] = [(align, name_fasta2)]
            list_score.append(score)
        else:
            dico_all[score].append((align, name_fasta2))
            list_score.append(score)
    list_score = reversed(sorted(list_score))
    make_output(name_query, list_score, dico_all, seq1, dico_psipred)

if __name__ == "__main__":
    FASTA_FILE = sys.argv[1]
    run_alignment_and_make_output(FASTA_FILE)





    #score_between_2_pssm(PSSM1,PSSM2)
    #mblos_path=CWD+"/data/BLOSUM62.txt"


    #sm = ScoringMatrix(seq1, seq2)
    #fill_matrix(sm,pssm1,pssm2)
    #align,score=get_alignments(sm,pssm1,pssm2)
    #make_output(name_fasta1,name_fasta2,align,score,seq1,seq2)
    #terminal_output.print_matrix(sm)
    #terminal_output.print_alignments(get_alignments(sm,pssm1,pssm2))
