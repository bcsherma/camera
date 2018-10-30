#!/usr/bin/env python
"""
checkPeaks.py

Comparing the peaks in 300ms vs 50ms experiments
"""

import argparse
import itertools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":

    # get command line arguments
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("hmqc")
    parser.add_argument("fifty_ms")
    parser.add_argument("three_hundred_ms")
    args = parser.parse_args()

    # get the csv files
    hmqc = pd.read_csv(args.hmqc)
    short = pd.read_csv(args.fifty_ms)
    long = pd.read_csv(args.three_hundred_ms)

    # x and y axes
    X = []
    Y = []

    # x and y geminal axes
    X_gem = []
    Y_gem = []

    # iterate over the rows of the 50 ms file
    for idx, entry in short.iterrows():

        # get coordinates
        c1 = entry["c1"]
        c2 = entry["c2"]
        h2 = entry["h2"]

        # ignore diagonals
        if abs(c1 - c2) < 0.15:
            continue

        # match receiver
        receiver = hmqc[(abs(hmqc["carbon"] - c2) < 0.2)
                        & (abs(hmqc["hydrogen"] - h2) < 0.02)]

        # match sender
        sender = hmqc[abs(hmqc["carbon"] - c1) < 0.2]

        if receiver.empty or sender.empty:
            continue

        # go over pairs of assignments
        for i, j in itertools.product(receiver["assignment"],
                                      sender["assignment"]):
            # check for geminality
            if i[:-1] == j[:-1]:
                geminal = True
                break

        else:  # not a geminal
            geminal = False

        # find the same NOEs in the 300ms file
        match = long[(abs(long["c1"] - c1) < 0.1)
                     & (abs(long["c2"] - c2) < 0.1)
                     & (abs(long["h2"] - h2) < 0.01)]

        if not match.empty:  # if there were matches, plot intensity

            # if np.mean(match["intensity"]) > 7*(10**7):
            #     print(c1, c2, h2)

            if geminal:
                X_gem.append(entry["intensity"])
                Y_gem.append(np.mean(match["intensity"]))

            else:  # otherwise ignore this
                X.append(entry["intensity"])
                Y.append(np.mean(match["intensity"]))

    # get upper and lower limits for axes
    upper = max(X + Y + X_gem + Y_gem)*(1.05)
    lower = min(X + Y + X_gem + Y_gem)*(0.9)

    plt.xlim(lower, upper)
    plt.ylim(lower, upper)
    plt.plot(np.linspace(lower, upper, 10), np.linspace(lower, upper, 10),
             "r--", label="diagonal")

    # make scatter plot
    plt.scatter(X, Y, label="non-geminal")
    plt.scatter(X_gem, Y_gem, c="g", label="geminal")
    plt.xlabel("50ms CCH intensity")
    plt.ylabel("300ms CCH intensity")
    plt.title("Intensity of NOEs in the 50 ms experiment vs the 300ms "
              "experiment")
    plt.legend()

    plt.show()
