#!/usr/bin/env python3
"""
meth2magma.py: converts methylphetamine spreadsheet into MAGMA input files

author: Ben Sherman
email: bcsherma@ucsc.edu
"""

import argparse
import pandas


def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__.strip())
    parser.add_argument("csv", type=str, help="csv file to be converted")
    parser.add_argument("output", type=str, help="name of output list")
    return parser.parse_args()


def main():
    """Main method for script"""
    args = get_args()
    csv = pandas.read_csv(args.csv)

    try:
        assert "cluster" in csv and "noe" in csv
    except AssertionError:
        print("ERROR! cluster and noe are required columns in the input csv")
        exit(1)

    noes = {(row["noe"][:row["noe"].rfind(".")],
             row["cluster"][:row["cluster"].rfind(".")])
            for _, row in csv.iterrows()}

    noes = {(a, b) for a, b in noes if a != b}  # filter out diagonals

    """
    Write the pairs to the output file
    """
    with open(args.output, "w") as outfile:
        for a, b in noes:
            outfile.write(f"{a} {b}\n")


if __name__ == "__main__":
    main()
