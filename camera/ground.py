"""
ground.py

Library of functions for checking whether or not ground truth is respected by
a given set of signatures, noe network, and structure.
"""

import itertools
import networkx as nx


def check_network(network, structure):
    """
    Given an NOE network, a collection of signatures, and a structure, verify
    that a maximum cardinality matching of the network does not necessarily
    violate the known assignments to the signatures
    """

    # Restrict the network to its active edges

    network = network.active_graph()

    # Iterate over connected components of this graph

    for component in nx.connected_component_subgraphs(network):
        
        min_radius = check_component(component, structure)
    
        if 10 < min_radius < float("inf"):

            print(f"warning: component cannot satisfy ground truth until "
                  f"{min_radius:.3f} angstroms")

            print("\t", end="")
            for c in component.nodes():
                print(c, end=" ")
            print("\n")


def check_component(component, structure):
    """
    Check whether any maximum cardinality matching of the given connected
    component of the symmetrization graph can respect ground truth
    """

    # Set the minimum observed distance of this component to be infinity

    minimum = float("inf")

    # Get the size of the maximum cardinality matching of this component

    max_matching_size = len(nx.max_weight_matching(component, 
                                                   maxcardinality=True)) 

    # Iterate over all subsets of edges of the component that are the same
    # size as the maximum cardinality matching

    for matching in itertools.combinations(component.edges(), 
                                           max_matching_size):

        # Check that this is in fact a matching

        if nx.is_matching(component, matching):

            # Use helper to compute the minimum length edge of this matching
            
            minimum = min(minimum, check_matching(matching, structure))

    # Return the minimum observed distance for this component

    return minimum
            

def check_matching(matching, structure):
    """
    Check whether a given matching, i.e. set of edges of the network, can be
    activated while respecting ground truth
    """

    return min([check_edge(i, j, structure) for i, j in matching])


def check_edge(alpha, beta, structure):
    """
    Check whether a given edge of the symmetrization graph violates ground
    truth
    """

    # Set the minimum distance to be infinity

    minimum = float("inf")
    
    # Iterate over all possible clusterings of the edge

    for alpha_cluster, beta_cluster in itertools.product(alpha.clusters,
                                                         beta.clusters):

        # Iterate over all possible assignments of these clusters
        
        for alpha_asg, beta_asg in itertools.product(alpha_cluster.asg,
                                                     beta_cluster.asg):

            # If these assignments are not the same, check the distance

            if alpha_asg != beta_asg:
                minimum = min(minimum, 
                              structure[alpha_asg][beta_asg]["distances"][0])

    # Return the minimum

    return minimum
     
