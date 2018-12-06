#!/usr/bin/env python3
"""
nmr2list.py

Converts nmr restraint file from RCSB into a list of pairs of methyl
identifiers in the conventional camera format, e.g. I1, L2.1
"""

import re
import argparse
import nmrstarlib

ATOM_EXPR = re.compile(r"M[GD][12]?")

COLORS = ["ALA", "ILE", "LEU", "VAL"]

SEQ1 = "Gen_dist_constraint.Seq_ID_1"
COMP1 = "Gen_dist_constraint.Comp_ID_1"
ATOM1 = "Gen_dist_constraint.Atom_ID_1"

SEQ2 = "Gen_dist_constraint.Seq_ID_2"
COMP2 = "Gen_dist_constraint.Comp_ID_2"
ATOM2 = "Gen_dist_constraint.Atom_ID_2"


def get_args():
    """
    Parse command line arguments
    """

    parser = argparse.ArgumentParser(description=__doc__.strip())
    parser.add_argument("restraints", help="nmr restraints file")
    parser.add_argument("-o", "--output", help="output file")

    return parser.parse_args()


def main():
    """
    Main method for the script
    """

    args = get_args()
    file = list(nmrstarlib.read_files(args.restraints))[0]
    restraints = file["save_CNS/XPLOR_distance_constraints_2"]["loop_1"]

    for r in restraints[1]:

        id1 = r[SEQ1]
        res1 = r[COMP1]
        atom1 = r[ATOM1]

        id2 = r[SEQ2]
        res2 = r[COMP2]
        atom2 = r[ATOM2]

        if res1 in COLORS and res2 in COLORS:
            if (ATOM_EXPR.match(atom1) is not None
                    and ATOM_EXPR.match(atom2) is not None):

                if res1 == "ILE":
                    label1 = f"I{id1}"
                else:
                    label1 = f"{res1[0]}{id1}.{atom1[-1]}"

                if res2 == "ILE":
                    label2 = f"I{id2}"
                else:
                    label2 = f"{res2[0]}{id2}.{atom2[-1]}"

                print(label1, label2)


if __name__ == "__main__":
    main()
