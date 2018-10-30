"""
symmetrize.py

Using the clustering CSP to determine if there are edges of the
symmetrization graph whose activation causes the formula to become
unsatisfiable. Such edges can be safely removed, which may cause more
connected components of size 2 to appear
"""

import tqdm
import halo
import itertools
import networkx as nx
from . import sat
from . import ground
from . import params

SPIN_STR = "checking a component of size {}..."


def reduce_symmetrization_graph(network, signatures, structure):
    """
    Given a network, signatures, and a structure, determine edges of the
    symmetrization graph to delete on the basis of them making the problem
    unsatisfiable.
    """

    # Run this outer loop until we make a pass over every connected component
    # of the symmterization graph and find that all edges in each component
    # are part of a maximum cardinality matching in which the problem is
    # satisfiable

    # Save how many edges are in the active part of the sym graph
    active_edges = len(network.active_graph().edges())

    iteration_no = 1
    while True:

        # Print a message about this iteration

        print(f"Beginning iteration number {iteration_no} over "
              f"components of the symmetrization graph\n")

        # Check for ground truth being respected

        ground.check_network(network, structure)

        # Iterate over living but inactive connected components of the
        # symmetrization graph.

        inactive_graph = network.inactive_graph()

        components = nx.connected_component_subgraphs(inactive_graph)
        components = sorted(components, key=lambda x: x.number_of_nodes())
        for component in components:

            # Call helper to determine which edges of this component are in
            # no maximum cardinality matchings that preserve SAT

            unseen = test_component(component, network, signatures, structure)

            # If there are edges which can deleted, do so and check to see
            # if this creates new components which are small enough to use as
            # distance constraints

            if unseen:

                old_active_edges = active_edges

                for i, j in unseen:
                    network.kill(i, j)

                # Reset the activity level and see if we got any new edges
                network.set_activity_level(params.MAX_COMP_SIZE)

                active_edges = len(network.active_graph().edges())

                # If we have new edges, restart the iteration

                if active_edges > old_active_edges:
                    break

        # If we were unable to make any improvements then stop
        else:
            print()
            break

        # Increase the iteration counter
        iteration_no += 1
        print()

    # Iterate over complex active components and kill edges which cannot be
    # activated

    ground.check_network(network, structure)
    clean_components(network, signatures, structure)


def clean_components(network, signatures, structure):
    """
    Iterate over active complex components of the symmetrization graph and
    kill all edges whose activation causes the problem to become UNSAT
    """

    # Print message

    print("Iterating over active complex components to simplify graph\n")

    # Construct satisfiability formula

    formula = sat.ClusteringCSP(signatures, network, structure)

    # Get the active components of the network

    active = network.active_graph()

    # Iterate over edges of the graph

    for i, j in tqdm.tqdm(active.edges()):

        # If these both have degree 1, then ignore

        if active.degree(i) == active.degree(j) == 1:
            continue

        # Force the edge to be active

        formula.add_aux_clause([formula.activation_variables[i][j]])

        # Run solver, and if the result in UNSAT kill i, j

        result = formula.solve()
        formula.flush()

        if not result:
            network.kill(i, j)

    # Print a newline
    print()


def test_component(component, network, signatures, structure):
    """
    Iterate over the maximum cardinality matchings of the given component. For
    each matching, test whether its activation causes satisfiability to be broken
    """

    text = SPIN_STR.format(component.number_of_nodes())

    spinner = halo.Halo(text=text, spinner='dots')
    spinner.start()

    unseen = set(component.edges())

    # Iterate over the max cardinality matchings of the graph

    for matching in max_matchings(component):

        # Activate the edges of the matching and test for satisfiability

        for i, j in matching:
            network.activate(i, j)

        if checksat(network, signatures, structure):

            for m in matching:
                if m in unseen:
                    unseen.remove(m)

        for i, j in matching:
            network.deactivate(i, j)

    if unseen:
        spinner.succeed(f"{text} removed {len(unseen)} edges")

    else:
        spinner.fail(f"{text} could not remove any edges")

    # Return the set of edges which were never observed

    return unseen


def checksat(network, signatures, structure):
    """
    Check whether the given network, signatures, and structure give rise to a
    satisfiable instance of the clutsering CSP
    """

    # Construct SAT formula

    formula = sat.ClusteringCSP(signatures, network, structure)

    # Run the solver

    result = formula.solve()

    if not result:
        return False

    else:
        return True


def max_matchings(graph):
    """
    Return an itertor over the maximum cardinality matchings of the graph
    """

    matching_size = len(nx.max_weight_matching(graph, maxcardinality=True))

    for matching in itertools.combinations(graph.edges(), matching_size):

        if nx.is_matching(graph, matching):

            yield matching
