#!/usr/bin/env python
"""
convert a 2D sparky list into a CSV file that contains the same information
"""

import argparse
import re
import pandas as pd


def get_args():
    """Get command line arguments"""
    # create parser
    parser = argparse.ArgumentParser()
    # add arguments to the parser
    parser.add_argument("list", help="2d sparky list")
    parser.add_argument("output", help="output CSV path")
    return parser.parse_args()


def parse_file(fname):
    """Parse sparky list into a series of triplets"""
    matcher = re.compile(r"^([AILVMT])([0-9]+)[CQH][BGED]([12]?)")

    # read in the file
    with open(fname) as file:
        lines = file.readlines()

    # initialize empty list of peaks
    peaks = []

    # iterate over input lines
    for idx, line in enumerate(lines):

        try:  # try to parse, print issue if it doesn't work

            label, carbon, hydrogen, *rest = line.split()  # split the line

            # get assignment information from the label
            color, id, order = matcher.match(label).groups()
            lv = color in ["L", "V"]
            label = f"{color}{id}.{order}" if lv else f"{color}{id}"

            # add of
            if lv:
                geminal = f"{color}{id}.{3 - int(order)}"
            else:
                geminal = ""

            asg = label + " " + geminal

            peaks.append([label, asg, geminal, color, float(carbon),
                          float(hydrogen)])

        # print where any issues happened
        except (AttributeError, ValueError):
            print(f"Bad peak on line {idx + 1}")

    # return peak list
    return peaks


def main():
    """Main method for the script"""
    args = get_args()
    peaks = parse_file(args.list)
    csv = pd.DataFrame(peaks, columns=["label", "assignment", "geminal",
                                       "color", "carbon", "hydrogen"])
    csv.to_csv(args.output, index=False)


# run main method for the script
if __name__ == "__main__":
    main()
