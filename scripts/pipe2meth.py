#!/usr/bin/env python
"""
Convert nmrpipe peak table into a csv with only the relevant rows
"""

import argparse
import itertools
import pandas as pd
import networkx as nx


class Peak:
    """
    Class for storing peaks so that we can eliminate duplicates using a
    connected component algorithm
    """

    def __init__(self, dict):
        """
        Construct from dictionary
        """

        # determine if this is an HCH or CCH peak
        self.type = "cch" if "c1" in dict else "hch"

        # set c1 and h1 depending on cch or hch
        self.c1 = float(dict["c1"]) if self.type == "cch" else 0.
        self.h1 = float(dict["h1"]) if self.type == "hch" else 0.

        # set ch and h2 coordinates
        self.c2 = float(dict["c2"])
        self.h2 = float(dict["h2"])

        # set intensity
        self.intensity = float(dict["intensity"])

    def connected(self, other, ctol, htol):
        """
        Determine if these are close enough in space to be connected
        """

        assert self.type == other.type  # only compare NOEs of same type

        if self.type == "cch":
            return (abs(self.c1 - other.c1) < ctol
                    and abs(self.c2 - other.c2) < ctol
                    and abs(self.h2 - other.h2) < htol)
        else:
            return (abs(self.h1 - other.h1) < htol
                    and abs(self.c2 - other.c2) < ctol
                    and abs(self.h2 - other.h2) < htol)


def component_to_dict(component, name):
    """Return a dictionary representation of single cluster of noes"""

    size = component.number_of_nodes()  # save size of component

    # initialize attributes
    c1 = 0.
    c2 = 0.
    h1 = 0.
    h2 = 0.
    intensity = 0.

    for node in component.nodes():
        c1 += node.c1/size
        c2 += node.c2/size
        h1 += node.h1/size
        h2 += node.h2/size
        intensity += node.intensity/size

    # return a blank dictionary
    return {"label": name, "c1": c1, "c2": c2, "h1": h1, "h2": h2,
            "intensity": intensity}


def get_args():
    """
    Construct argument parser and parse command line arguments
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("peaklist", type=str, help="pipe peaklist")
    parser.add_argument("output", type=str, help="output file")
    parser.add_argument("--ctol", type=float, default=0.1,
                        help="carbon tolerance from duplicate removal")
    parser.add_argument("--htol", type=float, default=0.01,
                        help="hydrogen tolerance from duplicate removal")
    return parser.parse_args()


def parse_pipe(fname, order=["h2", "c2", "c1"]):
    """Parse the important parts from a pipe file"""

    w1, w2, w3 = order  # get the order of axes

    with open(fname) as file:  # read in the pipe file
        lines = file.readlines()

    # find the line where the peaks declarations begin
    for idx, line in enumerate(lines):
        if line.startswith("VARS"):
            lines = lines[idx:]
            break

    # get variable names and entries
    vars = lines[0].split()[1:]
    entries = [e.split() for e in lines[2:]]
    new_entries = []

    # iterate over the entries of this table
    for idx, entry in enumerate(entries):

        new_entry = {"label": f"peak{idx}"}

        # iterate over attributes of the peak
        for v, e in zip(vars, entry):

            if v == "X_PPM":
                new_entry[w1] = e

            elif v == "Y_PPM":
                new_entry[w2] = e

            elif v == "Z_PPM":
                new_entry[w3] = e

            elif v == "HEIGHT":
                new_entry["intensity"] = e

        # add this new entry
        if len(new_entry.keys()) == 5:
            if float(new_entry["intensity"]) > 0:
                new_entries.append(Peak(new_entry))

    # return new peak list
    return new_entries


def main():
    """Main method for the script"""
    args = get_args()  # get command line arguments
    peaks = parse_pipe(args.peaklist)
    print(f"collected {len(peaks)} from table")

    final_list = []

    # construct a graph to eliminate duplicates
    graph = nx.Graph()
    graph.add_nodes_from(peaks)
    for i, j in itertools.combinations(peaks, 2):
        if i.connected(j, args.ctol, args.htol):
            graph.add_edge(i, j)

    # iterate over connected components of this graph
    for idx, cc in enumerate(nx.connected_component_subgraphs(graph)):
        final_list.append(component_to_dict(cc, f"peak{idx}"))

    csv = pd.DataFrame(final_list)  # construct dataframe
    csv.to_csv(args.output, index=False,
               columns=["label", "c1", "c2", "h2", "intensity"])


if __name__ == "__main__":
    main()
