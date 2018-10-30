#!/usr/bin/env python3
"""
sparky2meth.py:

Converts an assigned 3D sparky peak list into a methylphetamine formatted
.csv file with all columns filled out
"""

import argparse
import re
import pandas as pd


def get_methyl_ids(label):
    """
    Given a label from the 3D list, determine which methyls participate
    in the NOE and whether or not they comprise a geminal
    """
    matcher = re.compile(r"^([AILVMT])([1-9][0-9]*)[CQ][BDG]([12]?)"
                         r"-([AILVMT])?([1-9][0-9]*)?[CQ][BDG]([12]?)")

    c1, g1, o1, c2, g2, o2 = matcher.match(label).groups()  # run regex

    assert (c2 and g2) or (not c2 and not g2)

    if not c2:  # if no color, this is a geminal or diagonal
        c2 = c1
        g2 = g1

    # set labels according to color
    l1 = f"{c1}{g1}.{o1}" if c1 in ["L", "V"] else f"{c1}{g1}"
    l2 = f"{c2}{g2}.{o2}" if c2 in ["L", "V"] else f"{c2}{g2}"

    return l1, l2


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

    with open(args.sparkylist) as lines:

        lines = lines.readlines()
        for idx, line in enumerate(lines):
            idx = idx + 1

            try:
                label, c1, c2, h2, intensity, *rest = line.split()
                sender, receiver = get_methyl_ids(label)

            except ValueError:
                print(f"invalid NOE definition on line {idx}")
                continue

            except AttributeError:
                print(f"invalid NOE label on line {idx}")
                continue

            dic = {"cluster": receiver, "assignment": receiver, "noe": sender,
                   "label": label, "c1": c1, "c2": c2, "h2": h2,
                   "intensity": intensity, "color": receiver[0]}

            entries.append(dic)

    for idx, dic in enumerate(entries):

        """
        Iterate over NOEs to identify reciprocals
        """

        if "reciprocal" in dic:
            continue

        a = dic["noe"]

        b = dic["cluster"]

        for jdx, dic2 in enumerate(entries[idx+1:]):
            if dic2["cluster"] == a:
                if dic2["noe"] == b:

                    """
                    Point these NOEs to each other
                    """
                    entries[idx]["reciprocal"] = dic2["label"]
                    entries[idx + jdx + 1]["reciprocal"] = dic["label"]

                    """
                    If these are geminals, identify them as such in the
                    geminal column
                    """
                    if dic["cluster"][:-1] == dic2["cluster"][:-1]:
                        if dic["color"] in ["L", "V"]:
                            entries[idx]["geminal"] = dic2["label"]
                            entries[idx + jdx + 1]["geminal"] = dic["label"]

    csv = pd.DataFrame(entries)
    order = ["label", "assignment", "color", "cluster", "noe", "reciprocal",
             "geminal", "c1", "c2", "h2", "intensity"]

    csv.to_csv(args.output, columns=order, index=False)


if __name__ == "__main__":
    main()
