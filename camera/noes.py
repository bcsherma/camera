"""
noes.py

This file contains an object oriented representation of an NOE, either in
3D CCH/HCH form or 4D form. These are the observed interactions between the
2D peaks, or signatures, which are implemented in hmqc.py
"""

import pandas
import itertools
from . import params


class Noe:
    """
    Object oriented representation of an observed NOE, containing information
    about the position of the NOE and potentially information about the
    possible identities of the participants and reciprocals of the NOE
    """

    def __init__(self, source):
        """
        Construct an NOE object from a pandas series which contains all of 
        the relevant information about this NOE
        """

        # Convert the source object from a pandas series into a dictionary,
        # deleting all null entries

        source = source.dropna()
        source = source.to_dict()

        # Get necessary features from this NOE

        self.label = source.get("label", "")
        self.c2 = source["c2"]
        self.h2 = source["h2"]
        self.h1 = source.get("h1", None)
        self.c1 = source.get("c1", None)
        self.intensity = source.get("intensity", 0.)
        self.cluster_str = source.get("cluster", "").split()
        self.reciprocal_str = source.get("reciprocals", "").split()

        # Determine the type of this NOE

        if self.c1 and self.h1:
            self.type = "4D"

        elif self.c1:
            self.type = "CCH"

        else:
            self.type = "HCH"
        
        
        # Check for diagonals 

        if self.type == "CCH":
            if abs(self.c1 - self.c2) < 0.1:
                raise Warning
        
        elif self.type == "HCH":
            if abs(self.h1 - self.h2) < 0.01:
                raise Warning

        elif self.type == "4D":
            if abs(self.h1 - self.h2) < 0.01 and abs(self.c1 - self.c2) < 0.1:
                raise Warning

        # Initialize fields that will be filled in later

        self.clusters = []
        self.reciprocals = []

    def __repr__(self):
        """
        Return a text representation of this NOE for printing
        """

        return f"noe-{self.label}" 

    def __eq__(self, other):
        """
        Determine if this NOE and another are equal by comparing their labels
        """

        if not isinstance(other, Noe):
            return False

        return self.label == other.label

    def __hash__(self):
        """
        How to perform hashing on an NOE object
        """

        return hash(self.label)

    def symmetric(self, other):
        """
        Determine if this NOE and another are symmetric according to the 
        symmetrization tolerances
        """
        
        # If these are NOEs from different experiments, return False
        if self.type != other.type:
            return False

        elif self.type == "CCH":
            return (abs(self.c1 - other.c2) < params.SYM_CTOL 
                    and abs(self.c2 - other.c1) < params.SYM_CTOL)
        
        elif self.type == "HCH":
            return (abs(self.h1 - other.h2) < params.SYM_HTOL 
                    and abs(self.h2 - other.h1) < params.SYM_HTOL)
        
        # Check for 4D NOE Symmetry, which relies on all coordinates
        else:
            return (abs(self.h1 - other.h2) < params.SYM_HTOL  
                    and abs(self.h2 - other.h1) < params.SYM_HTOL
                    and abs(self.c1 - other.c2) < params.SYM_CTOL 
                    and abs(self.c2 - other.c1) < params.SYM_CTOL)

    def set_clusters(self, signatures):
        """
        Given a list of signatures, set this Noes clusters field to be all
        signatures which are within the clustering tolerance of itself
        """
        
        # Reset my list of clusters to be empty
        self.clusters = {}
        
        def clusterable(sig):
            return (abs(sig.carbon - self.c2) < params.CLS_CTOL
                    and abs(sig.hydrogen - self.h2) < params.CLS_HTOL)

        # Filter down to signatures within range
        signatures = filter(clusterable, signatures) 
        self.clusters = list(signatures)

    def to_dict(self):
        """
        Return a dictionary representation of this Noe
        """

        dictionary = {"label": self.label,
                      "intensity": self.intensity, 
                      "c1": self.c1, "c2": self.c2, "h2": self.h2,
                      "reciprocals": " ".join([
                          r.label for r in self.reciprocals]),
                      "clusters": " ".join([c.label for c in self.clusters])}
        
        if self.type == "4D":
            dictionary["h1"] = self.h1

        return dictionary

def set_clusters(noes, signatures):
    """
    Use the set clusters method of each instance of the Noe class
    in the given noe list.
    """

    # Iterate over the set of Noes and set the clusters of the Noe if it has
    # not already been manually clustered

    for noe in noes:

        if not noe.cluster_str:
            noe.set_clusters(signatures)

        else:
            noe.clusters = {s for s in signatures 
                            if s.label in noe.cluster_str}
        

def set_reciprocals(noes):
    """
    Given a set of NOEs, determine which have been set as each others
    reciprocals by force
    """
    
    pairs = itertools.combinations(noes, 2)
    pairs = [(i, j) for i, j in pairs 
             if i.label in j.reciprocal_str 
             or j.label in i.reciprocal_str]

    # Iterate over all pairs of known reciprocal and set them as such in the
    # reciprocal field

    for i, j in pairs:
        i.reciprocals.append(j)
        j.reciprocals.append(i)


def parse_noe_file(filename):
    """
    Given the name of an NOE CSV file, create a least of NOE objects
    and return them
    """

    noes = []
    csv = pandas.read_csv(filename)
    no_diagonals = 0 

    for idx, row in csv.iterrows():

        try:
            noes.append(Noe(row))
        
        except Warning:
            no_diagonals += 1

        except:
            print(f"Invalid NOE definition on line {idx + 1}"
                  f" of {filename}")

    # Return the list of Noe objects that we have aggregated
    print(f"Read {len(noes)} Noes from {filename} (excluding {no_diagonals} "
          f"diagonals)\n")
    return noes

