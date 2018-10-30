#!/usr/bin/env python
"""
peakmovement.py

Check how much the peaks are moved by centering by comaring the positions
of peaks in two different spreadsheets
"""

import pandas
import argparse


def get_args():
    """Parse command line aruments"""
    parser = argparse.ArgumentParser(description=__doc__.strip())
    # add arguments
    parser.add_argument("csv1", help="first csv file")
    parser.add_argument("csv2", help="second csv file")
    return parser.parse_args()  # return argument namespace


def main():
    """Main method for the script"""
    args = get_args()  # get command line arguments
    csv1 = pandas.read_csv(args.csv1)
    csv2 = pandas.read_csv(args.csv2)

    # iterate over the rows of the first csv file
    for idx, series in csv1.iterrows():

        # get all entries with the same label
        matches = csv2[csv2["label"] == series["label"]]

        # if unique, we will print the movement
        if len(matches) == 1:

            # get the unique match
            match = matches.iloc[0]

            # get the movement
            delta_c1 = abs(series["c1"] - match["c1"])
            delta_c2 = abs(series["c2"] - match["c2"])
            delta_h2 = abs(series["h2"] - match["h2"])
            delta_intensity = abs(series["intensity"] - match["intensity"])

            if delta_c1 < 0.15 and delta_c2 < 0.15 and delta_h2 < 0.02:
                continue

            # # print the difference
            # print(f"{series['label']:>20} {delta_c1:5.3f} {delta_c2:5.3f} "
            #       f"{delta_h2:5.3f} {delta_intensity:5.2E}")
            # print(f"{series['label']:>20} {series['c1']:5.3f} "
            #       f"{series['c2']:5.3f} {series['h2']:5.3f} "
            #       f"{series['intensity']:5.2E}")
            # print(f"{series['label']:>20} {match['c1']:5.3f} "
            #       f"{match['c2']:5.3f} {match['h2']:5.3f} "
            #       f"{match['intensity']:5.2E}")
            # print()
            print(series["label"])

if __name__ == "__main__":
    main()
