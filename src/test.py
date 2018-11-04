# This file is part of sequence-aligner.
# Copyright (C) 2014 Christopher Kyle Horton <chorton@ltu.edu>

# sequence-aligner is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# sequence-aligner is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with sequence-aligner. If not, see <http://www.gnu.org/licenses/>.


# MCS 5603 Intro to Bioinformatics, Fall 2014
# Christopher Kyle Horton (000516274), chorton@ltu.edu
# Last modified: 11/6/2014
import os
from scoring_matrix import ScoringMatrix
import numpy as np


gap_score = -1
terminal_gap_score = 0

def score_between_2_pssm(v1,v2):
    """
    this function goal is to calculate score bewteen 2 vectors

    """
    score=0
    cwd = os.getcwd()
    mblos_path=cwd+"/data/BLOSUM62.txt"
    m_blos=np.loadtxt(mblos_path)

    for i in range(0,len(v1)) :
        for j in range(0,len(v2)):
            score=score+v1[i]*v2[i]*m_blos[i][i]
    return score

def initialize_edges(sm, alignment_is_global=False):
    '''Sets up the top and left edges of the provided ScoringMatrix. This is
    the first step in the dynamic programming algorithm.
    Performs a semi-global alignment by default unless alignment_is_global is
    specified to be True.'''
    sm.set_score(0, 0, 0)
    gap = gap_score if alignment_is_global else terminal_gap_score
    for i in range(1, sm.get_rows()):
        sm.set_score(i, 0, sm.get_score(i - 1, 0) + gap)
        sm.add_up_backlink(i, 0)
    for i in range(1, sm.get_columns()):
        sm.set_score(0, i, sm.get_score(0, i - 1) + gap)
        sm.add_left_backlink(0, i)

def fill_matrix(sm,pssm1,pssm2, alignment_is_global=False):
    '''Uses dynamic programming to fill out a provided ScoringMatrix after
    the edges have been initialized.
    Performs a semi-global alignment by default unless alignment_is_global is
    specified to be True.
    Based on pseudocode provided on page 54 of our textbook.'''
    initialize_edges(sm, alignment_is_global)
    #print(sm.get_rows())
    #print(sm.get_columns())
    #print(pssm1.shape)
    #print(pssm2.shape)

    for i in range(1, sm.get_rows() +1):
        for j in range(1, sm.get_columns() +1):
            # Calculate scores
            score_diagonal = 0
            match=score_between_2_pssm(pssm1[i-1],pssm2[j-1])
            score_diagonal = sm.get_score(i - 1, j - 1) + match

            if (i == 0 or i == sm.get_rows() - 1) and not alignment_is_global:
                score_left = sm.get_score(i, j - 1) + terminal_gap_score
            else:
                score_left = sm.get_score(i, j - 1) + gap_score
            if (j == 0 or j == sm.get_columns() - 1) and not alignment_is_global:
                score_up = sm.get_score(i - 1, j) + terminal_gap_score
            else:
                score_up = sm.get_score(i - 1, j) + gap_score
            max_score = max(score_diagonal, score_left, score_up)
            sm.set_score(i, j, max_score)
            # Establish backlink(s)
            if max_score == score_diagonal:
                sm.add_diagonal_backlink(i, j)
            if max_score == score_left:
                sm.add_left_backlink(i, j)
            if max_score == score_up:
                sm.add_up_backlink(i, j)

def get_alignments(sm,pssm1,pssm2, alignment_is_global=False):
    '''Returns a list of the alignments generated from the scoring matrix.
    Performs a semi-global alignment by default unless alignment_is_global is
    specified to be True.
    Port of code from global-grid2.rb.'''
    fill_matrix(sm,pssm1,pssm2, alignment_is_global)
    # Find the optimal alignment paths starting from the lower right corner.
    seq = (sm.get_top_sequence(), sm.get_left_sequence())
    last_col = len(seq[0])
    last_row = len(seq[1])
    print(last_col)
    print(last_row)
    todo_list = [[last_row,last_col,"",""]] # Entry (row,col,string0,string1)
    print(todo_list)
    print(sm.get_backlinks(10, 17))
    done_list = [] # Entry (str0,str1)
    backlink_used = [[False for x in range(len(seq[0]) + 1)]
                      for x in range(len(seq[1]) + 1)]
    score=0
    while todo_list:
        row, col, str0, str1 = todo_list.pop()
        backlinks = sm.get_backlinks(row, col)
        if True in backlinks.values(): # If some back-link exists.
            score=score+sm.get_score(row,col)
            backlink_used[row][col] = True # Mark linked cells as used as we go.
            if backlinks["diagonal"]:
                todo_list.append([row - 1, col - 1, seq[0][col - 1] + str0,
                                 seq[1][row - 1] + str1])
            if backlinks["up"]:
                todo_list.append([row - 1, col, '_' + str0,
                                 seq[1][row - 1] + str1])
            if backlinks["left"]:
                todo_list.append([row, col - 1, seq[0][col - 1] + str0,
                                 '_' + str1])
        else:
            done_list.append([str0,str1])
    # Clean up unused backlinks
    for row in range(len(seq[1]) + 1):
        for col in range(len(seq[0]) + 1):
            if not backlink_used[row][col]:
                sm.remove_backlinks(row, col)

    print(score)
    return done_list

def read_pssm(path) :
    ##1ere ligne : nom du fasta
    ##2eme ligne : sequence
    ##a partir de la 3eme : tableau
    with open(path,"r") as fillin :
        name_fasta = fillin.readline()
        seq=fillin.readline().replace(" ","").replace("\n","")
        lines = (line for line in fillin)
        FH = np.loadtxt(lines)
    return(name_fasta,seq,FH)

if __name__ == "__main__":
    cwd = os.getcwd()
    name_fasta1,seq1,pssm1=read_pssm(cwd+"/data/query.aamtx")
    name_fasta2,seq2,pssm2=read_pssm(cwd+"/data/template.aamtx")
    #score_between_2_pssm(PSSM1,PSSM2)
    mblos_path=cwd+"/data/BLOSUM62.txt"

    n=pssm1.shape[0]
    m=pssm2.shape[0]
    seq_number1=""
    seq_number2=""
    for i in range(1,n) : seq_number1=seq_number1+str(i)
    for j in range(1,m) : seq_number2=seq_number2+str(j)

    import terminal_output

    sm = ScoringMatrix(seq1, seq2)
    fill_matrix(sm,pssm1,pssm2)
    #terminal_output.print_matrix(sm)
    terminal_output.print_alignments(get_alignments(sm,pssm1,pssm2))
