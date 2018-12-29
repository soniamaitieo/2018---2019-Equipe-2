"""
.. module:: Atom
   :synopsis: This module implements the Atom class.
"""


class Atom:
    """
    .. class:: Atom

      This class groups informations about an atom.

    Attributes:
        name (str): Name of the atom
        coords (Numpy array): 3D coordinates of the atom
    """

    def __init__(self, name):
        self.name = name
        self.coords = None

    def __str__(self):
        if self.coords is None:
            return "<" + self.name + " | " + "empty coordinates>"
        return "<" + self.name + " | " + str(self.coords) + ">"

    def set_coords(self, coords):
        """
            Sets the coordinates of the atom.

            Args:
                coords (Numpy array): x,y,z coordinates.
        """
        self.coords = coords
