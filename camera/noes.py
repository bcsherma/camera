"""
noes.py

This file contains an object oriented representation of an NOE, either in
3D CCH/HCH form or 4D form. These are the observed interactions between the
2D peaks, or signatures, which are implemented in hmqc.py
"""

import pandas


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
        self.reciprocal_str = source.get("reciprocal", "").split()

        # Determine the type of this NOE
        if self.c1 and self.h1:
            self.type = "4D"

        elif self.c1:
            self.type = "CCH"

        else:
            self.type = "HCH"
        
        # Initialize fields that will be filled in later
        self.clusters = {}
        self.reciprocal = []

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
            return (abs(self.c1 - other.c2) < 0.15 
                    and abs(self.c2 - other.c1) < 0.15)
        
        elif self.type == "HCH":
            return (abs(self.h1 - other.h2) < 0.01 
                    and abs(self.h2 - other.h1) < 0.01)
        
        # Check for 4D NOE Symmetry, which relies on all coordinates
        else:
            return (abs(self.h1 - other.h2) < 0.01 
                    and abs(self.h2 - other.h1) < 0.01
                    and abs(self.c1 - other.c2) < 0.15 
                    and abs(self.c2 - other.c1) < 0.15)

    def set_clusters(self, signatures):
        """
        Given a list of signatures, set this Noes clusters field to be all
        signatures which are within the clustering tolerance of itself
        """
        
        # Reset my list of clusters to be empty
        self.clusters = {}
        
        def clusterable(sig):
            return (abs(sig.carbon - self.c2) < 0.15 
                    and abs(sig.hydrogen - self.h2) < 0.02)

        # Filter down to signatures within range
        signatures = filter(clusterable, signatures) 
        self.clusters = set(signatures)



def set_clusters(noes, signatures):
    """
    Use the set clusters method of each instance of the Noe class
    in the given noe list.
    """

    # Iterate over the set of Noes and set the clusters of the Noe if it has
    # not already been manually clustered

    for noe in noes:
        if not noe.clusters:
            noe.set_clusters(signatures)


def parse_noe_file(filename):
    """
    Given the name of an NOE CSV file, create a least of NOE objects
    and return them
    """

    print("Reading NOEs from", filename)

    noes = []
    csv = pandas.read_csv(filename)
    
    for idx, row in csv.iterrows():
        try:
            noes.append(Noe(row))
        except:
            print(f"Invalid NOE definition on line {idx + 1}"
                  f" of {filename}")

    # Return the list of Noe objects that we have aggregated
    print()
    return noes
