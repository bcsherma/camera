#!/usr/bin/env python3
"""
Convert a 3D CCH peak file into a sparky list
"""

import pandas
import argparse


def get_args():
    """
    Parse command line arguments
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", help="3D CCH file (.csv)")
    parser.add_argument("output", help="output list (sparky format)")
    return parser.parse_args()


def main():
    """
    Main method for the script
    """

    # Parse command line arguments

    args = get_args()

    # Read in the CSV

    csv = pandas.read_csv(args.csv)

    # Write out each entry of the CSV to the list

    with open(args.output, "w") as outf:

        for _, peak in csv.iterrows():

            outf.write(f"?-?-?  {peak['c1']:.3f}   {peak['c2']:.3f}   "
                       f"{peak['h2']:.3f}\n")


if __name__ == "__main__":
    main()
