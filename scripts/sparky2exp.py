#!/usr/bin/env python3
"""
spark2exp.py
Converts a 2d sparky peaklist, e.g.

L5CD1-QD1     24.093      1.267
L5CD2-QD2     25.804      1.051
V9CG1-QG1     19.842      0.885
V9CG2-QG2     22.582      0.994
A11CB-QB      23.202      1.607
V12CG1-QG1    20.641      0.857
V12CG2-QG2    20.977      0.821

To .exp format, e.g.

9 ILE CD1 10.565
9 ILE HD1 0.983
13 ILE CD1 10.759
13 ILE HD1 0.869
20 ILE CD1 13.911
20 ILE HD1 0.625
55 ILE CD1 14.038
55 ILE HD1 0.939
66 ILE CD1 14.194
66 ILE HD1 0.224
77 ILE CD1 14.555
77 ILE HD1 0.803
117 ILE CD1 15.598
117 ILE HD1 0.668

author: Ben Sherman
email: bcsherma@ucsc.edu
"""


import re
import argparse


def main():
    """Main method for script"""

    """
    Construct argument parser and parse command line arguments
    """
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("peaklist", type=str, help="input sparky 2d peaklist")
    parser.add_argument("output", type=str, help="output file")
    args = parser.parse_args()

    """
    Regular expression for extracting relevant information from the peak labels
    """
    peakexpr = re.compile(r"([AILVMT])([1-9][0-9]*)"
                          r"([CQ][BDG][12]?)-([CQ][BDG][12]?)")

    """
    Dictionary mapping one letter color identifiers to three letter color
    identifiers
    """
    colormap = {"A": "ALA", "I": "ILE", "L": "LEU", "V": "VAL",
                "M": "MET", "T": "THR"}

    """
    Read in lines from input file and iterate over each one
    """
    outlines = []
    with open(args.peaklist) as lines:

        lines = lines.readlines()

        for idx, line in enumerate(lines):

            idx = idx + 1  # extract line number

            try:
                label, shift1, shift2, *rest = line.split()
                color, id, atom1, atom2 = peakexpr.match(label).groups()
                atom1 = atom1.replace("Q", "H")
                atom2 = atom2.replace("Q", "H")
                outlines.append(f"{id} {colormap[color]} {atom1} {shift1}")
                outlines.append(f"{id} {colormap[color]} {atom2} {shift2}")

            except ValueError:
                print("invalid peak definition on line", idx)

            except AttributeError:
                print("invalid peak definition on line", idx)

    """
    Write the output lines to a file
    """
    with open(args.output, "w") as outfile:
        outfile.write("\n".join(outlines) + "\n")


if __name__ == "__main__":
    main()
