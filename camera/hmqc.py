"""
hmqc.py

This file contains an object oriented representation of 2D NMR peaks to be
used in the methyl assignment problem.
"""

import pandas
import itertools


class Signature:
    """
    Class representation of a 2D NMR peak
    """

    def __init__(self, source):
        """
        Construct self from source, which should be a pandas series from a
        CSV file
        """

        # Convert the source object from a pandas series into a dictionary,
        # deleting all null entries
        source = source.dropna()
        source = source.to_dict()

        # Extract all relevant fields from the dictionary, or use default
        # value if that fails
        self.label = source["label"]
        self.carbon = source["carbon"]
        self.hydrogen = source["hydrogen"]
        self.color = source.get("color", "").split()
        self.asg_str = source.get("assignment", "").split()
        self.option_str = source.get("options", "").split()
        self.geminal_str = source.get("geminal", "")
        
        # Initialize fields that will be set later
        self.asg = []
        self.options = []
        self.geminal = []

    def is_geminal(self, other):
        """
        Determine if this Signature and another form a geminal pair based on
        their geminal field
        """

        return self.geminal == other

    def __eq__(self, other):
        """
        Determine if this Signature object and another are equal
        """

        return self.label == other.label

    def __hash__(self):
        """
        Hash a Signature object
        """

        return hash(self.label)

    def __repr__(self):
        """
        How to print out a Signature object
        """

        return f"signature-{self.label}"


def parse_hmqc_file(filename):
    """
    Given the name of an HMQC CSV file, create a least of Signature objects
    and return them
    """

    print("Reading HMQC peaks from", filename)

    signatures = []
    csv = pandas.read_csv(filename)
    
    for idx, row in csv.iterrows():
        try:
            signatures.append(Signature(row))

        except:
            print(f"warning: invalid HMQC peak definition on line {idx + 1} "
                  f"of {filename}")
    
    # Determine which pairs of signatures are known to be geminal pairs
    for i, j in itertools.combinations(signatures, 2):
        if i.geminal_str == j.label or j.geminal_str == i.label:

            # Set these to be each others geminal pair
            print(f"Setting {i} and {j} as a geminal pair")
            i.geminal = j
            j.geminal = i

    print()
    return signatures

