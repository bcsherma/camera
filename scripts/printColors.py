#!/usr/bin/env python3
"""
reads in as spreadsheet in which the clustering has been performed with color
labels and prints the number of clusters of each color
"""

import argparse
import pandas


def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", type=str, help="csv file to be read")
    return parser.parse_args()


def main():
    """Main method for script"""
    args = get_args()
    csv = pandas.read_csv(args.file, na_filter=False)
    colors = set(filter(None, csv["color"].tolist()))  # filter out empty str
    numbers = {c: 0 for c in colors}

    clusters = set(filter(None, csv["cluster"].tolist()))
    for c in clusters:
        numbers[c[0]] += 1

    for key in sorted(numbers.keys()):
        print(key, numbers[key])


if __name__ == "__main__":
    main()
