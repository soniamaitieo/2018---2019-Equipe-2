"""
.. module:: Template
   :synopsis: This module implements the Template class.
"""

# Third-party modulesimport logging
import logging
import numpy as np


logging.basicConfig(filename="run_warnings.log", level=logging.WARNING)


class Template:
    """
    .. class:: Template

      This class groups informations about a template sequence/structure.

    Attributes:
        name (str): Name of the template
        residues (list of Residue object): Template's sequence of residues as list of Residues
                                           objects
        benchmark (int): Fold family type of the template: 3: Family, 2: Superfamily,
                         1: Fold and 0: None
                         It tells how similar the template is from the query structure.
                         This is necessary to be able to benchmark the new scoring functions.
        pdb (str): PDB filename of the template
    """

    def __init__(self, name, residues):
        self.name = name            # ex: Agglutinin
        self.residues = residues
        self.benchmark = "."
        self.pdb = None             # ex: 1jlxa1
        self.reindexed_pdb = None   # ex: 1jlxa1_reindexed
        self.modeller_pdb = None    # ex: 1jlxa1_mod
        self.first = None

    def display(self):
        """
            This function prints the residues as a sequence

            Returns:
                str: Amino acid sequence.
        """
        return "".join(str(res.name) for res in self.residues)

    def set_benchmark(self, fold_type):
        """
            Sets the benchmark attribute of the template.

            Args:
                fold_type (int): 3: Family, 2: Superfamily, 1: Fold and 0: None
        """
        self.benchmark = fold_type

    def set_pdb_name(self, metafold_dict):
        """
            Set the pdb file name of the template from the template's name.

            Args:
                metafold_dict: A dictionary with key = template name and value = pdb file
        """
        self.pdb = metafold_dict[self.name].split(".")[0]

    def parse_pdb(self, pdb_path):
        """
            Parse the pdb file and set the CA coordinates.

            Args:
                pdb_path: Path to the directory containing .atm file of the actual template.
        """
        count_res = 0
        nb_atoms = 0
        flag = True
        with open(pdb_path, 'r') as file:
            for line in file:
                line_type = line[0:6].strip()
                name_at = line[12:16].strip()
                if line_type == "ATOM" and (name_at == "N" or name_at == "CA" or name_at == "C"):
                    res_id = int(line[22:26].strip())
                    # In some PDBs the residues ids do not start at 1
                    # so we remember the first id number
                    if flag:
                        self.first = res_id
                        flag = False
                    x_coord = float(line[30:38].strip())
                    y_coord = float(line[38:46].strip())
                    z_coord = float(line[46:54].strip())
                    if count_res <= len(self.residues):
                        # Skip gaps in the template
                        while count_res < len(self.residues) and \
                              self.residues[count_res].name == "-":
                            count_res += 1
                        if count_res == len(self.residues):
                            break
                        if line_type == "ATOM" and name_at == "N":
                            self.residues[count_res].n_atom\
                                .set_coords(np.array([x_coord, y_coord, z_coord]))
                            nb_atoms += 1
                        elif line_type == "ATOM" and name_at == "CA":
                            self.residues[count_res].ca_atom\
                                .set_coords(np.array([x_coord, y_coord, z_coord]))
                            nb_atoms += 1
                        elif line_type == "ATOM" and name_at == "CB":
                            self.residues[count_res].cb_atom\
                                .set_coords(np.array([x_coord, y_coord, z_coord]))
                            nb_atoms += 1
                        elif line_type == "ATOM" and name_at == "C":
                            self.residues[count_res].c_atom\
                                .set_coords(np.array([x_coord, y_coord, z_coord]))
                            nb_atoms += 1
                        if nb_atoms == 3:
                            count_res += 1
                            nb_atoms = 0

    def reindex_pdb_by_index(self, start_index=1, pdb_txt=''):
        """
        Original code from: https://zhanglab.ccmb.med.umich.edu/reindex_pdb/reindex_pdb.py
        Reindex residue number of PDB format text

        Args:
            start_index (int): Index of first residue
            pdb_txt (str): Text of input PDB to be reindexed

        Returns:
            str: Text of the reindexed PDB
        """
        pdb_txt_reindex = ''
        current_old_index = ''  # residue number in origin PDB
        warn_chain_id = ''  # warning about new chain ID

        for line in pdb_txt.splitlines():
            if len(line) < 27 or (not line.startswith("ATOM  ")\
                    and not line.startswith("HETATM") and not line.startswith("TER")):
                pdb_txt_reindex += line+'\n'
                continue
            elif not line[16] in ['A', ' ']:  # alternative location identifier
                continue
            res_seq = line[22:27]  # residue sequence number
            current_chain_id = line[21]  # chain identifier

            if not current_old_index:  # first residue encountered
                current_old_index = res_seq  # residue number in origin PDB
                current_new_index = int(start_index)
                chain_id = current_chain_id
                res_seq_new = str(current_new_index)
                res_seq_new = ' '*(4-len(res_seq_new))+res_seq_new+' '
            elif current_chain_id != chain_id:
                if warn_chain_id != current_chain_id:
                    logging.warning("Warning! Discarding chain '%s' in template '%s'.atm\n"
                                    % (current_chain_id, self.pdb))
                    warn_chain_id = current_chain_id
                continue
            elif res_seq != current_old_index:
                current_new_index += 1
                current_old_index = res_seq
                res_seq_new = str(current_new_index)
                res_seq_new = ' '*(4-len(res_seq_new))+res_seq_new+' '
            pdb_txt_reindex += line[:16]+' '+line[17:22]+res_seq_new+line[27:]+'\n'
        return pdb_txt_reindex

    def reindex_pdb(self, start_index, pdb_path, clean=True):
        """
            Original code from: https://zhanglab.ccmb.med.umich.edu/reindex_pdb/reindex_pdb.py
            
            * Parse PDB file
            * Reindex it according to start index
            * Write the new reindexed PDB

            Args:
                start_index (int): Index of first residue
                pdb_path (str): Path to the directory containing the PDB to reindex
        """
        pdb_file = pdb_path + "/" + self.pdb + ".atm"
        f_in = open(pdb_file, 'rU')
        pdb_txt = ''
        for line in f_in.read().splitlines():
            if line.startswith("END"):
                if clean:
                    line = line.replace("ENDMDL", "END   ")
                pdb_txt += line+'\n'
                break
            if (line.startswith("ATOM  ") or line.startswith("TER") or (not clean and not line[:6]\
                in ["DBREF ", "SEQADV", "MODRES", "HELIX ", "SHEET ", "SSBOND", "SITE  "])):
                pdb_txt += line+'\n'
        f_in.close()

        pdb_txt_reindex = self.reindex_pdb_by_index(start_index, pdb_txt)
        # Write the new PDB
        self.reindexed_pdb = self.pdb + "_reindexed"
        # Set the new name of the template's PDB
        pdb_file = pdb_path + "/" + self.reindexed_pdb + ".atm"
        with open(pdb_file, "w") as f_out:
            f_out.write(pdb_txt_reindex)
