#!/usr/bin/env python3
"""
test_clustering_csp.py

Run tests on the clustering CSP
"""

import os
import unittest
import networkx as nx
import camera.sat as sat
import camera.hmqc as hmqc
import camera.noes as noes
import camera.ground as ground
import camera.params as params
import camera.network as network
import camera.structures as structures
import camera.symmetrize as symmetrize


class TestClusteringCSP(unittest.TestCase):
    """
    Run tests on the clustering CSP
    """

    params.RADIUS = 10
    params.SHORT_RADIUS = 8

    test_dir = os.path.dirname(__file__) + "/hnh/"

    signatures = hmqc.parse_hmqc_file(test_dir + "hmqc.csv")

    crosspeaks = noes.parse_noe_file(test_dir + "merged.csv")

    noes.set_clusters(crosspeaks, signatures)
    noes.set_reciprocals(crosspeaks)

    crosspeaks = [c for c in crosspeaks if c.clusters]

    symgraph = network.SymGraph(crosspeaks, True)
    symgraph.ignore_geminals(signatures)
    symgraph.histogram()

    structure = structures.load_structure(test_dir + "g100.json")
    hmqc.set_assignment(signatures, structure)
    symgraph.set_activity_level(3)

    formula = sat.ClusteringCSP(signatures, symgraph, structure)

    def test_sat(self):
        """
        Test that the formula we have constructed is in fact satisfiable
        """

        self.assertTrue(TestClusteringCSP.formula.solve())

    def test_injection(self):
        """
        Test that each signature receives only one assignment, and that
        each methyl receives at most one assignment
        """

        solution = TestClusteringCSP.formula.solve()
        vertex_asg = [(alpha, beta) for vtype, alpha, beta in solution
                      if vtype == sat.Formula.ASG_VAR]

        for s in TestClusteringCSP.signatures:
            self.assertEqual(len([1 for a, b in vertex_asg if a == s]), 1)

        for m in TestClusteringCSP.structure.nodes():
            self.assertTrue(len([1 for a, b in vertex_asg if b == m]) <= 1)

    def test_respect_matching(self):
        """
        Test that a maximum cardinality matching is activated
        """

        active_graph = TestClusteringCSP.symgraph.active_graph()
        solution = TestClusteringCSP.formula.solve()
        active_edges = [(alpha, beta) for vtype, alpha, beta in solution
                        if vtype == sat.Formula.ACT_VAR]

        # Iterate over the active complex components of the symgraph

        for cc in nx.connected_component_subgraphs(active_graph):

            edges = [(a, b) for a, b in active_edges if a in cc.nodes()
                     and b in cc.nodes()]

            self.assertTrue(nx.is_matching(cc, edges))

    def test_clustering(self):
        """
        Test that one clustering is being enforced
        """

        active_graph = TestClusteringCSP.symgraph.active_graph()
        solution = TestClusteringCSP.formula.solve()
        clustering = [(alpha, beta) for vtype, alpha, beta in solution
                      if vtype == sat.Formula.CST_VAR]

        for node in active_graph.nodes():
            if len(node.clusters) > 1:
                self.assertEqual(len([1 for a, b in clustering if a == node]),
                                 1)

    def test_respect_active_edges(self):
        """
        Test that all active edges are being respected
        """

        active_graph = TestClusteringCSP.symgraph.active_graph()
        solution = TestClusteringCSP.formula.solve()

        assignment = {alpha: beta for vtype, alpha, beta in solution
                      if vtype == sat.Formula.ASG_VAR}

        clustering = {alpha: beta for vtype, alpha, beta in solution
                      if vtype == sat.Formula.CST_VAR}

        edges = [(i, j) for i, j in active_graph.edges()
                 if active_graph.degree(i) == active_graph.degree(j) == 1]

        edges += [(alpha, beta) for vtype, alpha, beta in solution
                  if vtype == sat.Formula.ACT_VAR]

        for i, j in edges:

            if len(i.clusters) == 1:
                i_sig = list(i.clusters)[0]

            else:
                i_sig = clustering[i]

            if len(j.clusters) == 1:
                j_sig = list(j.clusters)[0]

            else:
                j_sig = clustering[j]

            i_met, j_met = assignment[i_sig], assignment[j_sig]

            dist = TestClusteringCSP.structure[i_met][j_met]["distances"][0]

            if i.short_range:
                self.assertTrue(dist < 8)

            else:
                self.assertTrue(dist < 10)


if __name__ == "__main__":
    unittest.main()
