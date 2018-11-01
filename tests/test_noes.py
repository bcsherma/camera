#!/usr/bin/env python3
"""
test_clustering_csp.py

Run tests on the clustering CSP
"""

import pandas
import unittest
import camera.noes as noes


class TestNoes(unittest.TestCase):
    """
    Run tests on the NOE class
    """

    def test_diagonal(self):
        """
        Test that we ignore diagonals
        """

        noe = pandas.Series({"label": "p1", "c1": 21.1, "c2": 21.14,
                             "h2": 0.1})

        with self.assertRaises(Warning):
            noes.Noe(noe)

    def test_symmetry(self):
        """
        Test that symmetry checking works correctly
        """

        n1 = pandas.Series({"label": "p1", "c1": 21.1, "c2": 18.6, "h2": 0.1})
        n2 = pandas.Series({"label": "p1", "c1": 18.7, "c2": 21.2, "h2": 0.1})

        n1, n2 = noes.Noe(n1), noes.Noe(n2)

        self.assertTrue(n1.symmetric(n2))

        n2.c2 = 21.4

        self.assertFalse(n1.symmetric(n2))

        n2.c2 = 21.2
        n2.short_range = True

        self.assertFalse(n1.symmetric(n2))


if __name__ == "__main__":
    unittest.main()
