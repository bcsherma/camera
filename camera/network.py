"""
network.py

This file contains representation and methods for an NOE network used to
contrain the assignment of 2D peaks to methlyls for the methyl assignment
problem
"""

import pandas
import shutil
import itertools
import networkx as nx


class SignatureGraph(nx.Graph):
    """
    This class is a representation of the graph which has signatures as the
    vertices and edges between signatures that are known to have an Noe. One of
    these graphs represents and complete clustering and symmetrization of the
    NOE data.
    """

    def __init__(self, signatures, symgraph, clustering):
        """
        Construct an instance of the Signature graph class. First, call the
        super class constructor, then take all connected components of size 2
        from the symgraph and place an edge between the signatures to which
        they are mapped via the clustering dictionary
        """

        # Call the super class constructor
        nx.Graph.__init__(self)

        # Make my set of nodes the given signatures
        self.add_nodes_from(signatures)

        # Iterate over the connected components of size 2 in the living
        # symmetrization graph, adding edges to the graph accordingly

        living_network = symgraph.living_graph()

        for i, j in living_network.edges():

            if living_network.degree(i) == living_network.degree(j) == 1:
                self.add_edge(clustering[i], clustering[j], geminal=False)

        # Iterate over geminal pairs of signatures and add geminal edges
        # between them

        for i, j in itertools.combinations(signatures, 2):
            if i.is_geminal(j):
                self.add_edge(i, j, geminal=True)


class SymGraph(nx.Graph):
    """
    This class is an object representation of the symmetrization graph, which
    captures the potential symmetry, i.e. reciprocity of NOEs in a set of
    observations
    """

    def __init__(self, source, find_symmetries):
        """
        Construct a SymGraph object by first calling the constructor for
        the parent class, nx.Graph and then adding nodes and edges to self
        """

        # Call the nx.Graph constructor to initialize parent class fields
        nx.Graph.__init__(self)

        # Save the type of this dataset
        self.type = list({s.type for s in source})[0]

        # Add all noes from source as vertices of the graph
        self.add_nodes_from(source)

        # If we already know the symmetry between the NOES, then they are
        # stored in the reciprocals field of each NOE. Simply add the edges
        # to this graph between each NOE and their known reciprocals

        for i in self.nodes():
            for j in i.reciprocals:
                self.add_edge(i, j, active=False, dead=False)

        # Add an edge between all pairs of vertices which are symmetric. The
        # definition of symmetry is given in the Noe class in noes.py

        if find_symmetries:
            for i, j in itertools.combinations(source, 2):

                if i.reciprocals or j.reciprocals:
                    continue

                if i.symmetric(j):
                    self.add_edge(i, j, active=False, dead=False)

    def activate(self, i, j):
        """
        Activate the edge between i and j in the graph. If no such edge
        exists, raise an exception
        """

        if self.has_edge(i, j):
            self[i][j]["active"] = True

        else:
            raise ValueError(f"No edge between {i} and {j}")

    def deactivate(self, i, j):
        """
        Deactivate the edge between i and j in the graph. If no such edge
        exists, raise an exception
        """

        if self.has_edge(i, j):
            self[i][j]["active"] = False

        else:
            raise ValueError(f"No edge between {i} and {j}")

    def kill(self, i, j):
        """
        Kill the edge between i and j in the graph. If no such edge
        exists, raise an exception
        """

        if self.has_edge(i, j):
            self[i][j]["dead"] = True

        else:
            raise ValueError(f"No edge between {i} and {j}")

    def living_graph(self):
        """
        Return a copy of the graph with all of the edges that are dead removed
        """

        # Create a copy of this graph
        copy = self.copy()

        # Remove all edges which are dead
        copy.remove_edges_from([(i, j) for i, j in self.edges()
                                if self[i][j]["dead"]])

        # Return the copy with the dead edges removed
        return copy

    def active_graph(self):
        """
        Return a copy of the graph in which only active edges remain
        """

        # Get a copy of this graph with only the living edges
        copy = self.living_graph()

        # Remove all edges which are inactive
        copy.remove_edges_from([(i, j) for i, j in self.edges()
                                if not self[i][j]["active"]])

        # Remove vertices of degree 0 from the graph
        copy.remove_nodes_from([n for n in copy.nodes() if not copy.degree(n)])

        # Return the copy with only active edges
        return copy

    def inactive_graph(self):
        """
        Return a copy of the graph in which only inactive edges remain
        """

        # Get a copy of this graph with only the living edges
        copy = self.living_graph()

        # Remove all edges which are active
        copy.remove_edges_from([(i, j) for i, j in self.edges()
                                if self[i][j]["active"]])

        # Remove vertices of degree 0 from the graph
        copy.remove_nodes_from([n for n in copy.nodes() if not copy.degree(n)])

        # Return the copy with only active edges
        return copy

    def set_activity_level(self, max_size):
        """
        Set the all edges in connected components of the living graph with
        vertices less than or equal to some number to be active.
        """

        # Get the subgraph with living components
        living_graph = self.living_graph()

        # Iterate over connected components of the living graph
        for component in nx.connected_component_subgraphs(living_graph):

            # If the component has few enough vertices, active all edges.
            # Otherwise, deactivate all edges.
            if component.number_of_nodes() <= max_size:
                for i, j in component.edges():
                    self[i][j]["active"] = True

            else:
                for i, j in component.edges():
                    self[i][j]["active"] = False

    def to_csv(self, outfile):
        """
        Output all nodes as rows in a CSV file, where active edges are
        recorded as elements in the reciprocal column of the CSV
        """

        # Get a copy of the living component of the graph

        living_graph = self.living_graph()

        # Iterate over the nodes and them to running dict list

        dictionaries = []
        for n in self.nodes():

            if not n.reciprocals:
                if living_graph.has_node(n):
                    n.reciprocals = list(living_graph.neighbors(n))

            dictionaries.append(n.to_dict())

        # Write out to CSV file

        csv = pandas.DataFrame(dictionaries)

        # If the data type is 4D, we add the h1 column

        if self.type == "4D":
            columns = ["label", "c1", "h1", "c2", "h2", "intensity",
                       "clusters", "reciprocals"]

        else:
            columns = ["label", "c1", "c2", "h2", "intensity",
                       "clusters", "reciprocals"]

        # Write out the CSV file

        csv.to_csv(outfile, columns=columns, index=False)

    def histogram(self):
        """
        Print out a histogram of the sizes of the connected components of the
        living portion of the graph
        """

        # Get the living sym graph and the connected component sizes

        living_graph = self.living_graph()
        sizes = [component.number_of_nodes() for component
                 in nx.connected_component_subgraphs(living_graph)]
        all_sizes = set(sizes)
        counts = {s: sizes.count(s) for s in all_sizes}

        # Get the maximum line length from the size of the current terminal
        # window

        max_length = shutil.get_terminal_size()[0] - 30

        # Get the largest count and scale down all counts to proportion

        max_count = max(sizes.count(s) for s in all_sizes)

        if max_count > max_length:

            counts = {s: max(1, int(counts[s]*max_length/max_count))
                      for s in counts}

        # Iterate over all_sizes

        for size in sorted(all_sizes):

            print(f"[no.comp.of.size={size:<2}]:{sizes.count(size):<3}",
                  end="|")
            print(counts[size]*"\u25a7")

        print()
