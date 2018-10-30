#!/usr/bin/env python
"""
aggregateMap.py:

Takes MAPXS trial files and aggregates a score
"""

import argparse


def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__.strip())
    parser.add_argument("trials", type=str, nargs="+", help="trial files")
    return parser.parse_args()


def main():
    """Main method for the script"""
    args = get_args()
    trials = args.trials

    assignments = {}
    for t in trials:
        with open(t) as infile:
            lines = infile.readlines()
            start = None

            for idx, line in enumerate(lines):
                if line.startswith("Methyl Assignment Predictions:"):
                    start = idx + 1
                    break

            if not start:
                print("SHIT")
                exit()

            lines = lines[start:]
            for l in lines[:-1]:
                id, col, atom, _, id1, col1, atom1, *rest = l.split()
                h_vertex = "-".join([col, id, atom])
                g_vertex = "-".join([col1, id1, atom1])

                if h_vertex not in assignments:
                    assignments[h_vertex] = set()

                assignments[h_vertex].add(g_vertex)

    for a in assignments:
        print(a, assignments[a])

    total = len(assignments)
    correct = 0
    for a in assignments:
        for m in assignments[a]:
            if a[:-1] == m[:-1]:
                correct += 1
                break

    print(f"{100*correct/total}% of lists contain true assignment")


if __name__ == "__main__":
    main()
