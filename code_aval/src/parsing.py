"""
.. module:: parsing
   :synopsis: This module implements all the functions to parse either input
                or additional necessary files.
"""

# Third-party modules
import re
import os
import subprocess
import pathlib
import numpy as np
from Bio.SeqUtils import seq1

# Local modules
from src.residue import Residue
from src.alignment import Alignment
from src.query import Query
from src.template import Template


def parse_metafold(metafold_file):
    """
        Extracts name of metafold as key and associated pdb as
        value in dictionary.

        Args:
            metafold_file (str): A METAFOLD.list file containing for a given line a
                           template name and a pdb file name associated.

        Returns:
            dict: A dictionary with key = template name and value = pdb file.
    """
    metafold_dict = {}  # Initalization of the dictionary
    with open(metafold_file, "r") as file:
        for line in file:
            metafold_dict[line.split()[0]] = line.split()[1]
    return metafold_dict


def parse_foldrec(foldrec_file, metafold_dict):
    """
        Extracts the score, the template name, the query and template sequences for
        the first n alignments from a file containing N profil-profil alignments and
        creates a list of Alignment objects. It also gives the template's pdb file name
        and gets all the coordinates of the CA atoms in the template's Residue list.

        Args:
            foldrec_file (str): The file containing N profil-profil alignments and their
                                corresponding scores.
            metafold_dict (dict): A dictionary with key = template name
                                        and value = pdb file.

        Returns:
            dict: A dictionary with key = template name and value = an Alignment object.
    """
    # Regex :
    num_reg = re.compile("^No\\s*([0-9]+)")
    template_name_reg = re.compile("Alignment :.*vs\\s+([A-Za-z0-9-_]+)")
    score_reg = re.compile("^Score :\\s+([-0-9\\.]+)")
    query_seq_reg = re.compile("^Query\\s*([0-9]+)\\s*([A-Z0-9-]+)\\s*([0-9]+)")
    template_seq_reg = re.compile("^Template\\s*[0-9]+\\s*([A-Z-]+)")
    empty_query_reg = re.compile("^Query\\s+\\d\\s-+\\s+\\d.*$")

    alignment_dict = {}
    count_templates = 0

    with open(foldrec_file, "r") as file:
        query_reg_count = 0
        template_reg_count = 0
        for line in file:
            # Search a regex for the current line :
            num_found = re.search(num_reg, line)
            template_name_found = re.search(template_name_reg, line)
            score_found = re.search(score_reg, line)
            query_seq_found = re.search(query_seq_reg, line)
            template_seq_found = re.search(template_seq_reg, line)
            empty_query_found = re.search(empty_query_reg, line)
            # A num is found :
            if num_found:
                num = int(num_found.group(1))
            # A template name is found :
            if template_name_found:
                # These templates have more than 1 chain, so we must skip them
                if template_name_found.group(1) in ["bromodomain", "rhv", "Peptidase_A6", "ins",
                                                    "Arg_repressor_C", "SAM_decarbox", "prc",
                                                    "Chorismate_mut"]:
                    # We skip the alignment
                    for _ in range(10):
                        next(file)
                template_name = template_name_found.group(1)
            # A score is found :
            if score_found:
                score = float(score_found.group(1))
            # Empty alignement (the query = gaps only)
            if empty_query_found and query_reg_count == 2:
                print("Skipping alignement "+str(num)+". The query is only composed of gaps.")
                # We skip the alignment
                for _ in range(5):
                    next(file)
            # A query sequence is found :
            if query_seq_found:
                if query_reg_count == 0:
                    query_first = int(query_seq_found.group(1))
                    query_seq = [Residue(name) for name in list(query_seq_found.group(2))]
                    query_last = int(query_seq_found.group(3))
                    query_reg_count += 1
                elif query_reg_count == 1:
                    for ind, sec_struct in enumerate(list(query_seq_found.group(2))):
                        query_seq[ind].secondary_struct = sec_struct
                    query_reg_count += 1
                elif query_reg_count == 2:
                    for ind, ss_conf in enumerate(list(query_seq_found.group(2))):
                        if ss_conf != "-":
                            ss_conf = int(ss_conf)
                        query_seq[ind].ss_confidence = ss_conf
                    query_reg_count = 0
            # A template sequence is founds :
            if template_seq_found:
                if template_reg_count == 0:
                    template_seq = [Residue(name) for name in list(template_seq_found.group(1))]
                    template_reg_count += 1
                elif template_reg_count == 1:
                    for ind, sec_struct in enumerate(list(template_seq_found.group(1))):
                        template_seq[ind].secondary_struct = sec_struct
                    template_reg_count = 0
                    # Add a new alignment object in the list :
                    ali = Alignment(num, score,
                                    Query(query_seq, query_first, query_last),
                                    Template(template_name, template_seq))
                    ali.template.set_pdb_name(metafold_dict)
                    ali.template.parse_pdb("data/pdb/"+ali.template.name+"/"+ali.template.pdb
                                           +".atm")
                    alignment_dict[template_name] = ali
                    count_templates += 1
    return alignment_dict


def parse_benchmark(query_name, benchmark_file, alignment_dict):
    """
        Extract the line of the benchmark file containing the benchmarks of
        the studied query. Each benchmark template have an assigned structure
        corresponding to the degree of similarity with the query
        ("Family", "Superfamily" or "Fold"). For a given benchmark template
        object, the structure is stored in the benchmark item.

        Args:
            query_name (str): Name of the query.
            benchmark_file (str): The path to the benchmark file.
            alignment_dict (dict): A dictionary with key = template name.
                                   and value = an Alignment object.
    """
    with open(benchmark_file, "r") as file:
        for line in file:
            # Only the line containing the query name is extracted
            if query_name in line:
                # List containing for each benchmark the template name and
                # the associated fold_type
                template_fold_type = line[0:-1].split(": ")[1].split(", ")
                for i in template_fold_type:
                    # The fold_type of the benchmark template is stored in
                    # the benchmark item
                    template = i.split(" ")[0]
                    fold_type = i.split(" ")[1]
                    if template in alignment_dict:
                        alignment_dict[template].template.set_benchmark(fold_type)
                # No need to continue because all benchmarks are stored
                break


def parse_dope(dope_file):
    """
        Extracts 30 dope energy values for the 20*20 residus-CA pairs and stores it
        in a dictionary with key = res_1-res_2 and value = an array of the 30 dope scores.

        Args:
            dope_file (str): The file dope.par containing energy values for each
                       amino acid pair.

        Returns:
            dict: A dictionary with key = res_1-res_2 and value = an array of
            the 30 dope energy values.
    """
    dope_dict = {}
    with open(dope_file, "r") as file:
        for line in file:
            # get the line with C-alpha for both amino acids
            if line[4:6] == "CA" and line[11:13] == "CA":
                res_1 = seq1(line[0:3])
                res_2 = seq1(line[7:10])
                dope_dict[res_1+res_2] = np.fromstring(line[14:-1], dtype=float, sep=" ")
    return dope_dict


def get_index_list(aln_file_clustal):
    """
        Get the index list of the query.

        Args:
            aln_file_clustal (str): Path to the multiple alignment file (clustal format).

        Returns:
            list: List of the index of the query with gaps.
    """
    # First line of this multiple alignment
    with open(aln_file_clustal, "r") as file:
        query_seq = file.readline()[:-1]
        # "non_gap" positions index are saved into a list
        index_list = [i for i, ele in enumerate(query_seq) if ele != "-"]
    return index_list


def parse_ccmpred_result(ccmpred_result, index_list):
    """
        Create a dictionary of the top couplings of the ccmpred result file.
        Usefull for the co-evolution score calculation.

        Args:
            ccmpred_result (str): Path to the output file generated by ccmpred.
            index_list (list): List of indexes corresponding to amino acids positions in aln file.            

        Returns:
            dict: A dictionary with key = ranking of coupling based on ss_confidence
            and value = index aa1, index aa2
    """
    # Number of top coupling expected
    ntops = int(len(index_list)/2)
    # Extract tops coupling
    top_couplings = subprocess.check_output(
        ["./bin/CCMpred/scripts/top_couplings.py -n {} {}".
         format(str(ntops), ccmpred_result)], shell=True).decode('utf-8').split('\n')[1:-1]
    # Create a dictionary of top couplings
    top_couplings_dict = {}
    for k, value in enumerate(top_couplings):
        values = value.split("\t")
        index_i = int(values[0])
        index_j = int(values[1])
        # Do not parse gaps associated with top couplings
        if (index_i in index_list) and (index_j in index_list):
            top_couplings_dict[k] = (index_list.index(index_i), index_list.index(index_j))
    return top_couplings_dict
