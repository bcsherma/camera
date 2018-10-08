"""
network.py

This file contains representation and methods for an NOE network used to
contrain the assignment of 2D peaks to methlyls for the methyl assignment
problem
"""

import itertools
import networkx as nx


class SymGraph(nx.Graph):
    """
    This class is an object representation of the symmetrization graph, which
    captures the potential symmetry, i.e. reciprocity of NOEs in a set of 
    observations
    """
    
    def __init__(self, source, symmetry_known):
        """
        Construct a SymGraph object by first calling the constructor for 
        the parent class, nx.Graph and then adding nodes and edges to self
        """

        # Call the nx.Graph constructor to initialize parent class fields
        nx.Graph.__init__(self)

        # Add all noes from source as vertices of the graph 
        self.add_nodes_from(source)
        
        # If we already know the symmetry between the NOES, then they are
        # stored in the reciprocals field of each NOE. Simply add the edges
        # to this graph between each NOE and their known reciprocals
        if symmetry_known:

            for i in self.nodes():
                
                for j in i.reciprocals:
                    self.add_edge(i, j, active=False, dead=False)

        else:
            # Add an edge between all pairs of vertices which are symmetric. The 
            # definition of symmetry is given in the Noe class in noes.py
            for i, j in itertools.combinations(source, 2):
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

        # Remove vertices of degree 0 from the graph
        copy.remove_nodes_from([n for n in copy.nodes() if not copy.degree(n)])

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
