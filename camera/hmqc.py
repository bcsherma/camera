"""
hmqc.py

This file contains an object oriented representation of 2D NMR peaks to be
used in the methyl assignment problem.
"""

import pandas
import shutil
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
        self.color = list(source.get("color", ""))
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

    def to_dict(self):
        """
        Convert self into dictionary
        """

        if self.asg:
            self.asg_str = " ".join([m.label for m in self.asg])

        if self.options:
            self.option_str = " ".join([m.label for m in self.options])

        return {"label": self.label,
                "color": "".join(sorted(self.color)),
                "assignment": " ".join(self.asg_str),
                "options": " ".join(self.option_str),
                "carbon": f"{self.carbon:.3f}",
                "hydrogen": f"{self.hydrogen:.3f}",
                "geminal": self.geminal.label if self.geminal else ""}

    def __eq__(self, other):
        """
        Determine if this Signature object and another are equal
        """

        if not isinstance(other, Signature):
            return False

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

    def nailed(self):
        """
        Return True if my options are a single methyl or a geminal pair
        """

        return len({m.seqid for m in self.options}) == 1


def set_assignment(signatures, structure):
    """
    Given a list of signatures and a structure, set the assignment and options
    fields of each signature so that each respectively contains the set of
    methyls compatible with the given assignment and options from the
    CSV file
    """

    # Get a set containing the nodes of the structure
    methyls = set(structure.nodes)

    # Iterate over the signatures in the list
    for sig in signatures:

        # Update the assignment and options lists to contain the methyls
        # with labels in the asg_str and options_str fields respectively
        sig.asg = {m for m in methyls if m.label in sig.asg_str}
        sig.options = {m for m in methyls if m.label in sig.option_str}


def to_csv(signatures, outfile):
    """
    Write out set of signatures to csv file
    """

    dictionaries = [s.to_dict() for s in signatures]
    csv = pandas.DataFrame(dictionaries)
    csv.to_csv(outfile, columns=["label", "color", "assignment", "options",
                                 "geminal", "carbon", "hydrogen"], index=False)


def nailed_histogram(signatures, support=None):
    """
    Print a histogram of support set sizes, either from given support
    dictionary or using the options field of the hmqc peaks
    """

    # If support is None, then set it as a mapping from each signature to
    # its options field

    if support is None:
        support = {s: s.options for s in signatures}

    # Get support set sizes

    sizes = [len(support[s]) for s in support if support[s]]
    all_sizes = set(sizes)

    # Get maximum histogram bar size
    max_length = shutil.get_terminal_size()[0] - 20

    # Iterate over all_sizes
    for size in sorted(all_sizes):

        print(f"[no.options={size:<2}]:{sizes.count(size):<3}", end="|")
        print(min(sizes.count(size), max_length) * "\u25a7")

    # Figure out how many nailed there are

    total = len(support)

    nailed = 0
    for sig, sup in support.items():

        if len(sup) == 1:
            nailed += 1

        elif len(sup) == 2:

            left, right = list(sup)
            if left.geminal(right):
                nailed += 1

    print(f"\n{100*nailed/total:.3f}% nailed\n")


def parse_hmqc_file(filename):
    """
    Given the name of an HMQC CSV file, create a least of Signature objects
    and return them
    """

    print("Reading HMQC peaks from", filename, "\n")

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

    # Print out a histogram of the number of types

    print()
    type = ["".join(sorted(s.color)) for s in signatures]
    all_types = set(type)

    for t in sorted(all_types):
        print(f"Read {type.count(t)} signatures of color {t}")

    print()
    return signatures
