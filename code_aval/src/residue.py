"""
.. module:: Residue
   :synopsis: This module implements the Residue class.
"""

# Third-party modules
import numpy as np

# Local modules
from src.atom import Atom


class Residue:
    """
    .. class:: Residue

      This class groups informations about a residue.

    Attributes:
        name (str): Name of the residue (1 letter code)
        ca_coords (Numpy array): 3D coordinates of the residue
    """

    def __init__(self, name):
        self.name = name
        self.ca_atom = Atom("CA")
        self.cb_atom = Atom("CB")
        self.c_atom = Atom("C")
        self.n_atom = Atom("N")
        self.secondary_struct = None
        self.ss_confidence = None

    def __str__(self):
        return str(self.name)

    def __repr__(self):
        return str(self.name)

    def calculate_distance(self, residue, carbon):
        """
            Calculate Euclidian distance between two residues with THE most efficient method, the
            Einstein summation convention with *numpy.einsum*.

            .. math::

               distance = \sqrt{(x_a-x_b)^2+(y_a-y_b)^2+(z_a-z_b)^2}

            Args:
                residue (object): An object of the Residue class.
                carbon (str): Type of the carbon.

            Returns:
                float: The calculated distance.
        """
        if carbon == "CA":
            a_min_b = self.ca_atom.coords - residue.ca_atom.coords
        elif carbon == "CB":
            a_min_b = self.cb_atom.coords - residue.cb_atom.coords
        # The Numpy's "einsum" function is the most efficient way to calculate a distance
        #  accordingly to this benchmarking: https://stackoverflow.com/a/47775357/6401758
        dist = np.sqrt(np.einsum('i,i->', a_min_b, a_min_b))
        return dist
