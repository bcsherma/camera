#!/usr/bin/env python3
"""
sparky2csv.py:

Converts an un-assigned 3D sparky peak list into a methylphetamine formatted
.csv file with the basic columns filled out
"""

import argparse
import pandas as pd


def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__.strip())
    parser.add_argument("sparkylist", type=str,
                        help="sparky file to be converted")
    parser.add_argument("output", type=str, help="name of output csv file")
    parser.add_argument("--hch", action="store_true",
                        help="hch instead of cch")
    return parser.parse_args()


def main():
    """Main method for the script"""
    args = get_args()

    entries = []

    noe_dim = "h1" if args.hch else "c1"  # save the name of the noe dimension

    with open(args.sparkylist) as lines:

        lines = lines.readlines()
        # lines = set(lines)  # remove duplicate lines
        peak = 1
        for idx, line in enumerate(lines):
            idx = idx + 1

            try:
                label, c1, c2, h2, intensity, *rest = line.split()

                c1 = float(c1)  # convert these to floats
                c2 = float(c2)
                h2 = float(h2)
                intensity = float(intensity)

                label = f"peak{peak}"
                peak += 1

            except ValueError:
                print(f"invalid NOE definition on line {idx}")
                continue

            dic = {"label": label, noe_dim: c1,
                   "c2": c2, "h2": h2, "intensity": intensity}

            entries.append(dic)

    # create dataframe and write out
    csv = pd.DataFrame(entries)
    order = ["label", noe_dim, "c2", "h2", "intensity"]
    csv.to_csv(args.output, columns=order, index=False)


if __name__ == "__main__":
    main()
