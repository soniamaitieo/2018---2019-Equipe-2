"""
.. module:: Query
   :synopsis: This module implements the Query class.
"""

from src.residue import Residue

class Query:
    """
    .. class:: Query

      This class groups informations about a query sequence.

    Attributes:
        residues (list of Residue objects): Template's sequence of residues as list of Residues
                                            objects
        first (int): First residue of the query sequence.
        last (int): Last residue of the query sequence.
    """

    def __init__(self, residues, first, last):
        self.residues = residues
        # This is a "temporary" list of residues for Modeller, which we will modify only
        # in order to make the PIR alignment (by inserting gaps).
        self.modeller_residues = residues
        self.first = first
        self.last = last

    def display(self, modeller=False):
        """
            This function prints the residues as a sequence

            Args:
                modeller (boolean): True/False

            Returns:
                str: Amino acid sequence.
        """
        if modeller:
            return "".join(str(res.name) for res in self.modeller_residues)
        return "".join(str(res.name) for res in self.residues)

    def get_size(self):
        """
            Get the size of the query sequence.

            Returns:
                int: The size of the query sequence.

        """
        return len(self.residues)

    def add_gaps_in_query(self, template):
        """
        This functions inserts the same gaps that occur in the template, in the query.

        Args:
            template (Template object): The current alignment's template
        """
        # Find the indexes of gaps in the template
        occ_ind = [i for i, res in enumerate(template.residues) if res.name == "-"]
        for i in occ_ind:
            # Insert the gaps in the query
            self.modeller_residues[i:i] = [Residue("-")]
