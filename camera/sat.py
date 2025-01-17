"""
sat.py

Library for creating, modifying and solving propositional satisfiability
formulae for methyl assignment.
"""

import os
import re
import uuid
import tqdm
import random
import itertools
import subprocess
import networkx as nx
from . import params


# Regular expression for getting assignments from cryptominisat
SOL_EXPR = re.compile(r"(-?[1-9][0-9]*)")


class Formula:
    """
    Basic formula class. Here the methods for solving a formulae are
    contained
    """

    # Enumerating variable types
    ASG_VAR = 1
    CST_VAR = 2
    ACT_VAR = 3
    CMD_VAR = 4
    EDG_VAR = 5

    def __init__(self):
        """
        Create a basic instance of the formula class
        """

        # Initialize the fields that all formulae have
        self.nvars = 0
        self.nclauses = 0
        self.base_clauses = []
        self.aux_clauses = []

        # Create variable meaning table
        #
        # Variable meaning is a table which maps boolean varibles to their
        # semantics. The semantics of a variable are reprsented a 3-tuple.
        # The first element of the tuple says whether the variable is an
        # assignment, clustering, activation, or commander variable. If the
        # variable is an assignment variable, the second element is the
        # relevant signature and the third is the relevenant methyl. If the
        # variable is a clustering variable, the second element is the
        # relevant NOE and the third is the relevant signature. If it is an
        # activation variable, then the 2nd a 3rd elements are the relevant
        # NOEs (in arbitrary order).

        self.variable_meaning = {}

    def add_clause(self, lits):
        """
        Add dysjunction of given literals to the formula
        """

        self.nclauses += 1
        self.base_clauses.append(lits)

    def add_aux_clause(self, lits):
        """
        Add an auxiliary, i.e. temporary clause to the fomrula
        """

        self.nclauses += 1
        self.aux_clauses.append(lits)

    def flush(self):
        """
        Remove all auxiliary clauses from the formula
        """

        self.nclauses = self.nclauses - len(self.aux_clauses)
        self.aux_clauses = []

    def next_variable(self):
        """
        Get the next available variable and increment the number of variables
        that have been assigned
        """

        self.nvars += 1
        return self.nvars

    def enumerate(self):
        """
        Enumerate the support sets for the signatures
        """

        # Print out a message

        print("Beginning the enumeration of support sets\n")

        # Localize the assignment variable table

        asgvar = self.assignment_variables

        # Initialize set of all vertices as those which have not had their
        # support sets fully enumerated

        unfinished = set(asgvar.keys())
        pbar = tqdm.tqdm(total=len(unfinished))

        # Initialize support sets to be empty

        support = {u: set() for u in unfinished}

        # Loop until we are done with the enumeration

        while True:

            # If unfinished is now empty, break

            if len(unfinished) == 0:

                # Close the progress bar and leave a newline
                pbar.close()
                print()

                # Return the computed support sets
                return support

            # Select an unfinished vertex at random

            focus = random.sample(unfinished, 1)[0]

            # Force this vertex to take an as yet unseen assignment

            for seen in support[focus]:
                self.add_aux_clause([-asgvar[focus][seen]])

            # Run the solver

            result = self.solve()
            self.flush()  # Delete aux clauses

            if result:

                # Update support sets

                for vtype, alpha, beta in result:
                    if vtype == Formula.ASG_VAR:
                        support[alpha].add(beta)

            else:

                # We are done with this vertex, so remove it from unfinished
                # and lock it to its known assignments
                pbar.update()
                unfinished.remove(focus)
                self.add_clause([asgvar[focus][s] for s in support[focus]])

    def inject_vertices(self, signatures, structure):
        """
        Create assignment variables for each signature-methyl pair such that
        the signature is not forbidden from being assigned to the methyl
        on the basis of support sets, forced assignment, or color. Constrain
        the assignment variables so that exactly one is true for any given
        signature or methyl, i.e. make the mapping in any satisfying
        assignment an injective function
        """

        # Get the set of methyls from the structure
        methyls = list(structure.nodes())

        # Iterate over the signatures
        for signature in signatures:

            # Localize the assignment variable table for this signature
            table = self.assignment_variables[signature]

            # Check for whether we are meant to filter by assignment and
            # support sets

            if params.FORCE_SV and signature.options:
                domain = signature.options

            elif params.FORCE_ASG and signature.asg:
                domain = signature.asg

            else:
                # Otherwise, just filter down by color compatiblity
                domain = filter(lambda m: m.color in signature.color, methyls)

            # Create a variable for each methyl in the domain
            for methyl in domain:

                # Get a variable for this methyl assignment
                var = self.next_variable()
                table[methyl] = var
                self.variable_meaning[var] = (Formula.ASG_VAR, signature,
                                              methyl)

            # Force exactly one of the assignment variables for this signature
            # to be True
            lits = [table[m] for m in table]
            self.at_most_one(lits)
            self.add_clause(lits)

        # Force each methyl to have no more than one signature assigned to it
        for m in methyls:

            # Get signatures that can be assigned to this methyl
            domain = [s for s in signatures if m in
                      self.assignment_variables[s]]

            # Force no more than one of the assignment variables for
            # assignments to m to be true in any satisfying assignment
            self.at_most_one([self.assignment_variables[d][m]
                              for d in domain])

    def at_most_one(self, lits):
        """
        Use the commander encoding to force that at most one of the given
        literals is True in any satisfying assignment to the formula
        """

        # If there are 3 or fewer literals in the given list, simply
        # use the naive implmentation of at_most_one
        if len(lits) < 4:
            self.naive_at_most_one(lits)
            return

        # Initialize an empty list of commander variables
        commanders = []

        # Partition the given literals into small sets. Create a commander
        # variable for each set and force the commander to be True IFF one
        # variable in its set is true
        idx = 0
        while idx < len(lits):

            # Get partition of variables and create commander for it
            group = lits[idx:idx+3]
            cmdr = self.next_variable()
            commanders.append(cmdr)

            # Set the meaning of the commander variable
            self.variable_meaning[cmdr] = (Formula.CMD_VAR, None, None)

            # Make the commander imply that one of the literals is true
            self.add_clause([-cmdr] + group)

            # Make the commander imply that every literal in its set if false
            for l in group:
                self.add_clause([cmdr, -l])

            # Make it so that at most one of the literals in the group is true
            self.at_most_one(group)

            idx += 3

        # Make it so that at most one of the commanders we created
        # is true in any satisfying assignment
        self.at_most_one(commanders)

    def naive_at_most_one(self, lits):
        """
        Use the naive clause construction to force no more than one of the
        given literals to be True in any satisfying assignment
        """

        for l1, l2 in itertools.combinations(lits, 2):
            self.add_clause([-l1, -l2])

    def to_string(self):
        """
        Return a DIMACS format string of this formula
        """

        # Get the header of the formula and the body separately
        header = f"p cnf {self.nvars} {self.nclauses}"
        body = "\n".join(map(clause_to_string,
                             self.base_clauses + self.aux_clauses))

        # Return the header and the body combined
        return header + "\n" + body

    def to_file(self, outfile):
        """
        Write formula out to output file
        """

        with open(outfile, "w") as outf:
            outf.write(self.to_string())

    def solve(self):
        """
        Run the solver and get back a solution
        """

        # Get string representation of self
        formula_str = self.to_string()

        # Run the solver as a subprocess
        process = subprocess.run(["cryptominisat5", "--verb=0"],
                                 input=bytes(formula_str, "utf-8"),
                                 stdout=subprocess.PIPE,
                                 timeout=15)

        assignments = get_assignments(process.stdout)

        return [self.variable_meaning[a] for a in assignments if a > 0]


def get_assignments(output):
    """
    Parse a complete variable assignment from the output of cryptominisat
    """

    # Use REGEX to find all integer strings in the output string and convert
    # to ints and then return.
    return set(map(int, SOL_EXPR.findall(str(output))))


def clause_to_string(clause):
    """
    Return a string representation of the given clause in DIMACS format
    """

    return " ".join(map(str, clause)) + " 0"


class ClusteringCSP(Formula):
    """
    This class extends formula to implement the clustering CSP, which takes
    as input a collection of signatures, an NOE network, and a structure, then
    constructs a CSP which constraint the assignment of signatures into the
    structure using the network to form ambiguous constraints.
    """

    def __init__(self, signatures, network, structure):
        """
        Construct a new clustering CSP given signatures, an NOE network, and
        a structure. First, call the Formula constructor to get define the
        basic methods which make use of the SAT solver.
        """

        # First, call the Formula constructor
        Formula.__init__(self)

        # Replace the network with the active components of the network
        network = network.active_graph()

        # Create variable tables
        #
        # self.assignment_variables is a dictionary where the keys are
        # signatures and the items are dictionaries. The item dictiories have
        # methyls as keys and boolean variable identifiers, i.e. integers as
        # items. self.assignment_variables[s][m] returns the boolean variable
        # representing the proposition "s is assigned to m".
        #
        # self.clustering_variables is the same idea but it is a mapping from
        # NOEs to dictioanries which map signatures to boolean variables.
        # self.clustering_variables[n][s] returns the boolean variable
        # representing the proposition "n is clustered to s".
        #
        # self.activation_variables maps noes to mappings from their neighbors
        # to boolean variables. self.activation_variables[n][n2] returns the
        # boolean variable representing the proposition that (n, n2) are truly
        # reciprocal NOEs.

        self.assignment_variables = {s: {} for s in signatures}
        self.clustering_variables = {n: {} for n in network.nodes()}
        self.activation_variables = {n: {} for n in network.nodes()}

        # Construct the formula by creating variables and clauses

        self.inject_vertices(signatures, structure)
        self.create_clustering_variables(network, signatures)
        self.create_activation_variables(network)
        self.respect_matching(network)
        self.distance_constraints(signatures, network, structure)
        self.geminal_constraints(signatures, structure)

    def enumerate_clusterings(self):
        """
        Enumerate the clusterings of the data
        """

        # Initialize an empty list of clusterings. Each element of the list
        # is a dictionary mapping NOEs to the signatures they are clustered to
        # in a solution. Each clustering is initialize with all uniquely
        # clusterable NOEs as keys mapping to their unique clusters

        clusterings = []

        while True:

            # Run the solver and get a solution back

            solution = self.solve()

            # If UNSAT, then flush aux clauses from the formula and return
            # all the clusterings we found so far

            if not solution:
                self.flush()
                return clusterings

            # Iterate over the clustering variables set to true by in the
            # discovered solution. Forbid this clustering from reoccuring and
            # add it to the list of found clusterings

            clause = []
            clustering = {}
            for node in self.clustering_variables.keys():
                if len(node.clusters) == 1:
                    clustering[node] = list(node.clusters)[0]

            for vtype, node, cluster in solution:
                if vtype == Formula.CST_VAR:
                    clustering[node] = cluster
                    clause.append(-self.clustering_variables[node][cluster])

            self.add_clause(clause)
            clusterings.append(clustering)

    def create_clustering_variables(self, network, signatures):
        """
        Create variables that represent the clustering of the NOES in the NOE
        network. Make it so that each NOE gets precisely one clustering.
        """

        # Iterate over the nodes of the NOE network
        for noe in network.nodes():

            # If this noe has only one possible clustering, do not create
            # a variable because its clustering does not vary

            if len(noe.clusters) == 1:
                continue

            # Localize the clustering variable table this NOE

            table = self.clustering_variables[noe]

            # Iterate over clusters this NOE can be clustered to

            for cluster in noe.clusters:

                # Create a variable that captures the clustering of this noe
                # to this cluster

                var = self.next_variable()
                self.variable_meaning[var] = (Formula.CST_VAR, noe, cluster)
                table[cluster] = var

            # Make it so that exactly one of the clustering variables is true

            lits = [table[c] for c in table]
            self.at_most_one(lits)
            self.add_clause(lits)

    def create_activation_variables(self, network):
        """
        Create variables corresponding to the activation of individual edges
        in the NOE network.
        """

        # Iterate over edges of the network

        for i, j in network.edges():

            # If this edge is its own connected component of the graph, it is
            # always active. Therefore we do not need to create a variable for
            # its activity.

            if network.degree(i) == 1 and network.degree(j) == 1:
                continue

            # Create a variable representing the activity of this edge

            var = self.next_variable()
            self.variable_meaning[var] = (Formula.ACT_VAR, i, j)
            self.activation_variables[i][j] = var
            self.activation_variables[j][i] = var

    def respect_matching(self, network):
        """
        Force a maximum cardinality matching of each connected component of
        the NOE network to be respected by each satisfying assignment to the
        formula.
        """

        # Iterate over the connected components of the network

        for component in nx.connected_component_subgraphs(network):

            # If the component has only two nodes it is respected by default

            if component.number_of_nodes() < 3:
                continue

            # Partition the nodes into biparite sets

            left, right = nx.bipartite.sets(component)

            # As a convention, the left side is the one with fewer vertices.
            # Swap left and right if right has fewer vertices.

            if len(left) > len(right):
                left, right = right, left

            # The way we encode the constraint that a maximum cardinality
            # matching of the symmetrization graph be respected is to force
            # one edge incident to each vertex on the left side of each
            # component be activated, while at most one edge incident to each
            # vertex of the right side be activated. This will activate a
            # maximum cardinality matching of the graph IFF the size of the
            # maximum cardinality matching of each component is the number of
            # vertices in the smaller biparite set of that component. Here we
            # run a check to verify that this is the case.

            mcm_size = len(nx.max_weight_matching(component,
                                                  maxcardinality=True))

            assert mcm_size == len(left)

            # Iterate over vertices on the smaller biparite set and force
            # exactly one of their incident edges to be activated

            for l in left:
                table = self.activation_variables[l]
                lits = [table[r] for r in table]
                self.at_most_one(lits)
                self.add_clause(lits)

            # Iterate over vertices on the larger bipartite set and force
            # at most one of their incident edges to be activated

            for r in right:
                table = self.activation_variables[r]
                lits = [table[l] for l in table]
                self.at_most_one(lits)

    def distance_constraints(self, signatures, network, structure):
        """
        Force all satisfying assignments to respect distance constraints
        given the activated edges and clustering.
        """

        # Iterate over edges of the network

        for i, j in network.edges():

            short = i.short_range

            # Iterate over clusterings of this edge

            for i_c, j_c in itertools.product(i.clusters, j.clusters):

                # If i_c and j_c are the same, then we can't possible respect
                # the constraint, so ignore

                if i_c == j_c:
                    continue

                # base_clause is the negation of the conjuction of the
                # variable indicating that this edge is active with the
                # variables saying that this is the clustering. If both
                # clusterings are unique and the edge is always active, then
                # the base clause is empty.

                base_clause = []

                if network.degree(i) > 1 or network.degree(j) > 1:
                    base_clause.append(-self.activation_variables[i][j])

                if len(i.clusters) > 1:
                    base_clause.append(-self.clustering_variables[i][i_c])

                if len(j.clusters) > 1:
                    base_clause.append(-self.clustering_variables[j][j_c])

                # Call a helper that forces an edge between i_c and j_c to
                # be respected, conditional on the base clause being
                # unsatisfied, i.e. conditional on (i,j) being active, i being
                # clustered to i_c and j being clustered to j_c

                self.respect_distance_constraint(i_c, j_c, structure,
                                                 base_clause, short)

    def respect_distance_constraint(self, alpha, beta, structure, base_clause,
                                    short):
        """
        Force an edge between alpha and beta (signatures) to be respected,
        conditional on the given base clause being unsatisfiable
        """

        # Get the alpha and beta assignment variable tables locally
        alpha_table = self.assignment_variables[alpha]
        beta_table = self.assignment_variables[beta]

        # Iterate over the domain of alpha

        for alpha_methyl in alpha_table:

            # Create a copy of the base_clause
            clause = base_clause.copy()

            # Append negation of alpha -> alpha_methyl variable to clause
            clause.append(-alpha_table[alpha_methyl])

            # Localize the neighborhood of alpha_methyl in the structure
            alpha_neighborhood = structure[alpha_methyl]

            # Iterate over the domain of beta
            for beta_methyl in beta_table:

                # If alpha_methyl and beta_methyl are the same, then this
                # can't be an edge of G and can't be used to satisfy the
                # edge between alpha and beta

                if alpha_methyl == beta_methyl:
                    continue

                # If alpha_methyl and beta_methyl are close in the structure,
                # then allow alpha -> alpha_methyl and beta -> beta_methyl
                # to satisfy the clause

                distance = alpha_neighborhood[beta_methyl]["distances"][0]

                if short:
                    if distance < params.SHORT_RADIUS:
                        clause.append(beta_table[beta_methyl])

                elif alpha_methyl.added or beta_methyl.added:
                    if distance < params.ADDED_RADIUS:
                        clause.append(beta_table[beta_methyl])

                elif distance < params.RADIUS:
                    clause.append(beta_table[beta_methyl])

            # Add this clause to the formula
            self.add_clause(clause)

    def geminal_constraints(self, signatures, structure):
        """
        Force all geminal pairs of signatures to be assigned to geminal pairs
        in the structure by every satisfying assignment to the formula
        """

        # Identify all geminal pairs of signatures

        geminals = {(i, i.geminal) for i in signatures if i.geminal}

        # Iterate over geminals pairs of signatures

        for i, j in geminals:

            # Iterate over the domain of i

            for i_methyl in self.assignment_variables[i]:

                clause = [-self.assignment_variables[i][i_methyl]]

                # Should i be assigned to i_methyl, the clause we add
                # can only be satisfied if j is assigned to the geminal pair
                # of i_methyl

                for j_methyl in self.assignment_variables[j]:
                    if i_methyl.geminal(j_methyl):
                        clause.append(self.assignment_variables[j][j_methyl])

                # Add this clause to the formula
                self.add_clause(clause)


class IsomorphismCSP(Formula):
    """
    CSP which maps one graph into another, without modelling ambiguity of
    clustering or symmetrization
    """

    def __init__(self, graph_h, structure, edge_variables=True):
        """
        Construct an instance of the subgraph isomorphism csp, first
        calling the super class constructor
        """

        # Call super class constructor to set up fundamental methods and fields

        Formula.__init__(self)

        # If edge_variables is set to False, we do not create variables for
        # the edge assignments and we do not allow for sampling.

        self.edge_vars = edge_variables

        # Set up variable tables.
        #
        # The assignment variables table maps the vertices of the graph H to
        # a dictionary which maps the methyls in the domain of a vertex to the
        # variable representing that assignment.
        #
        # The variable_cost maps variables to the cost of having that variable
        # set to True. This is used to store the cost of each edge assignment
        # so that we can perform weighted sampling.

        self.assignment_variables = {n: {} for n in graph_h.nodes()}
        self.variable_cost = {}

        # Set up variables and clauses

        self.inject_vertices(graph_h.nodes(), structure)
        self.distance_constraints(graph_h, structure)

    def marginalize(self, exponent, num_samples):
        """
        Get an empirical estimate of the marginal probability of each
        assignment using the SPUR sampler
        """

        # Get a list of samples

        samples = self.sample(exponent, num_samples)

        # Initialize a table that stores the marginal probability of every
        # assignment

        marginals = {n: {} for n in self.assignment_variables.keys()}

        # Iterate over the samples and update the marginals

        for sample in samples:

            for vertex, assignment in sample.items():

                if assignment.seqid in marginals[vertex]:
                    marginals[vertex][assignment.seqid] += 1/num_samples

                else:
                    marginals[vertex][assignment.seqid] = 1/num_samples

        # Return the marginals and the samples

        return marginals, samples

    def sample(self, exponent, num_samples):
        """
        Run the SPUR weighted sampling code to obtain the desired number of
        samples. The probability of each satisfying assignment is 1 over the
        product of the length of the used edges of G, raised to the given
        exponent
        """

        # Assert that we have created edge variables

        assert self.edge_vars

        # Generate a random filename for the formula and weightfile that will
        # be used for the sampler

        basename = "/tmp/" + str(uuid.uuid4())

        # Write out the weights file

        with open(basename + ".txt", "w") as weightfile:
            for variable, cost in self.variable_cost.items():
                weightfile.write(f"{variable} 1 {cost**-exponent}\n")

        # Write the cnf file out to a temporary file

        self.to_file(basename + ".cnf")

        # Run the sampler

        subprocess.run(f"spur -s {num_samples} -cnf {basename + '.cnf'} "
                       f"-w {basename + '.txt'} "
                       f"-out {basename + '.samples'}",
                       shell=True, stdout=open(os.devnull, "w"))

        # Read in the samples

        with open(basename + ".samples") as file:
            lines = file.readlines()
            lines = list(filter(lambda x: x.strip(), lines))
            lines = list(filter(lambda x: x[1] == ",", lines))
            samples = list(map(lambda x: x.strip()[2:], lines))
            assignments = self.parse_assignments(samples)

        # Delete all the temporary files

        os.remove(basename + ".cnf")
        os.remove(basename + ".txt")
        os.remove(basename + ".samples")

        # Return the assignments we sampled

        return assignments

    def parse_assignments(self, samples):
        """
        Parse assignments string output by SPUR sampler
        """

        # Iterate over the samples we collected and parse the assignments.
        # Each assignment is a mapping from signatures to methyls

        assignments = []
        for sample in samples:
            assignment = {}
            for idx, val in enumerate(sample):
                if int(val):
                    vtype, alpha, beta = self.variable_meaning[idx + 1]
                    if vtype == Formula.ASG_VAR:
                        assignment[alpha] = beta
            assignments.append(assignment)

        # Return the list of assignments

        return assignments

    def distance_constraints(self, graph_h, structure):
        """
        Force each edge of the graph_h to be assigned to an edge of the
        structure. Geminal edges of the graph H must be assigned to geminal
        edges of the structure.
        """

        # Localize the assignment variable table

        asgvar = self.assignment_variables

        # Iterate over the edges of the graph H

        for i, j in graph_h.edges():

            # Check whether these are a geminal
            h_geminal = i.is_geminal(j)

            # Iterate over the domain of the vertex i

            for i_met in asgvar[i].keys():

                # Initialize an empty clause which needs to be satisfied by
                # the assignment of j if i is assigned to i_met

                clause = [-asgvar[i][i_met]]

                # Iterate over domain of the vertex j

                for j_met in asgvar[j].keys():

                    # If i_met and j_met are the same methyl, ignore

                    if i_met == j_met:
                        continue

                    # Compare geminality of the edge of H to this pair in the
                    # structure

                    if h_geminal > i_met.geminal(j_met):
                        continue

                    # Get the distance between i_met and j_met over the
                    # structure

                    distance = structure[i_met][j_met]["distances"][0]

                    # If the distance between i_met and j_met is short enough,
                    # allow j -> j_met to satisfy this clause

                    within_range = False

                    if graph_h[i][j]["short"]:
                        if distance < params.SHORT_RADIUS:
                            within_range = True

                    elif i_met.added or j_met.added:
                        if distance < params.ADDED_RADIUS:
                            within_range = True

                    elif distance < params.RADIUS:
                        within_range = True

                    if within_range:

                        clause.append(asgvar[j][j_met])

                        if not self.edge_vars:
                            continue

                        # Create a variable which represents that (i,j) is
                        # mapped to i_map, j_map

                        variable = self.next_variable()
                        self.variable_cost[variable] = distance
                        self.variable_meaning[variable] = (Formula.EDG_VAR,
                                                           (i, j),
                                                           (i_met, j_met))

                        # Constrain the variable to be true if and only if the
                        # vertex assignment variables are true

                        self.add_clause([-variable, asgvar[i][i_met]])
                        self.add_clause([-variable, asgvar[j][j_met]])
                        self.add_clause([variable, -asgvar[i][i_met],
                                         -asgvar[j][j_met]])

                # Add the clause to the set of base clauses

                self.add_clause(clause)


def gibbs_reduce(graph_h, structure, exponent, num_samples):
    """
    Sample from the gibbs distribution and rerun enumeration, then rinse and
    repeat until we reach a fixed point
    """

    # Construct an isomorphism CSP from the given graph_h and structure

    formula = IsomorphismCSP(graph_h, structure)

    old_percent_nailed = 0.

    unnailed = {n for n in graph_h.nodes() if not n.nailed()}

    # Loop forever

    while True:

        # Initialize the percent nailed to be zero.

        percent_nailed = 0.

        # Run marginalization.

        marginals, samples = formula.marginalize(exponent, num_samples)
        correlate_assignments(formula, marginals, samples)

        for vertex, distribution in marginals.items():

            if vertex not in unnailed:
                continue

            highest = max(distribution.keys(), key=lambda x: distribution[x])

            if distribution[highest] > 0.9:

                # Increase the percent nailed for this vertex and nail it to
                # its assignments

                percent_nailed += 1/graph_h.number_of_nodes()
                options = [m for m in structure.nodes() if m.seqid == highest
                           and m in formula.assignment_variables[vertex]]
                formula.add_clause([formula.assignment_variables[vertex][m]
                                    for m in options])
                unnailed.remove(vertex)

                print(f"Nailing {vertex} to {options}")

        # If we did not increase the percent nailed from the previous iteration
        # then stop now

        if percent_nailed == old_percent_nailed:
            print()
            return formula.enumerate()

        # We increase the percent nailed, so save that progress

        old_percent_nailed = percent_nailed


def correlate_assignments(formula, marginals, samples):
    """
    Iterate over all pairs of vertex assignments. If the probability of a pair
    co-occuring deviates significantly from the product of the marginal
    probability of those assignments, then place a hard constraint forcing
    the assignments to be either logically equivalent or mutually exclusive.
    """

    # Localize the formula's assignment variable table.

    asgvar = formula.assignment_variables

    # Iterate over all pairs of signatures

    for sig1, sig2 in itertools.combinations(marginals.keys(), 2):

        # If these are a geminal pair, ignore.

        if sig1.is_geminal(sig2):
            continue

        # If either of these are nailed, ignore.

        if len(marginals[sig1]) == 1 or len(marginals[sig2]) == 1:
            continue

        # Iterate over all pairs of assignments between these two.

        for m1_id, m2_id in itertools.product(marginals[sig1].keys(),
                                              marginals[sig2].keys()):

            # If these are the same methyl, ignore.

            if m1_id == m2_id:
                continue

            # Get the product of the marginal probability of sig1 -> m1 and
            # sig2 -> m2.

            independent_prob = marginals[sig1][m1_id] * marginals[sig2][m2_id]

            # Count the fraction of samples in which these assignments co-occur

            real_prob = 0
            for sample in samples:

                if sample[sig1].seqid == m1_id and sample[sig2].seqid == m2_id:
                    real_prob += 1/len(samples)

            # Get the methyls with the given sequence IDs.

            m1s = [m for m in asgvar[sig1].keys() if m.seqid == m1_id]

            m2s = [m for m in asgvar[sig2].keys() if m.seqid == m2_id]

            if real_prob > 10*independent_prob:

                for m1 in m1s:
                    clause = [-1*asgvar[sig1][m1]]
                    for m2 in m2s:
                        clause.append(asgvar[sig2][m2])
                    formula.add_clause(clause)

                for m2 in m2s:
                    clause = [-1*asgvar[sig2][m2]]
                    for m1 in m1s:
                        clause.append(asgvar[sig1][m1])
                    formula.add_clause(clause)

            elif real_prob*10 < independent_prob:
                for m1, m2 in itertools.product(m1s, m2s):
                    formula.add_clause([-asgvar[sig1][m1], -asgvar[sig2][m2]])
