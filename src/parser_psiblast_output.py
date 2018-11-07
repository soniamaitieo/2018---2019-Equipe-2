#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import sys

def parse_psiblast_output(psiblastoutput , matchfile): 
    #Open psiblast file
    with open(psiblastoutput ,'r') as f:
        lines = f.readlines()     
    # Only take matchs from third round
    for num,line in enumerate(lines, 0):
        if line.startswith("Results from round 3"):
            index = num
            break
    #Create list of all match            
    match_query_list=[]
    for x in range( index , len(lines)) :
        if lines[x].startswith(">") :
            match_query_list.append(x)
    #Create dico with match-seq
    dico={}
    for value in range(0,len(match_query_list) - 1):
        v1 = match_query_list[value]
        v2= match_query_list[value + 1]
        L=[]
        all_scores=[]
        for x in range(v1,v2):
            if lines[x].startswith(" Score") :
                all_scores.append(x)
        if len(all_scores) > 1 :
            v0 = match_query_list[value]
            v1 = all_scores[0]
            v2 = all_scores[1]
            for x in range(v1,v2):
                if lines[x].startswith("Sbjct") :
                    L.append(lines[x].split()[2])
                    print(L)
            dico[lines[v0].split()[1]] = ''.join(L)
        else :
            for x in range(v1,v2):
                if lines[x].startswith("Sbjct") :
                    L.append(lines[x].split()[2])
            dico[lines[v1].split()[1]] = ''.join(L)
    with open( matchfile ,'w') as fout:
        for k, v in dico.items():
            fout.write(''.join([">" , k , "\n" , v , "\n"]))


if __name__ == "__main__":
    blastoutput = sys.argv[1]
    matchfile = sys.argv[2]
    parse_psiblast_output(blastoutput, matchfile)
