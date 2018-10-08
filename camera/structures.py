"""
structure.py

Module for reading pdb structures and storing information about the methyls
and the distances between all methyl pairs
"""

import json
import Bio.PDB
import itertools
import networkx as nx

# Create global pdb parsing object
PARSER = Bio.PDB.PDBParser(QUIET=True)


def load_structure(filename):
    """
    Load a graphical representation of the structure from a json file
    """

    # Read in json dictionary
    with open(filename) as infile:
        dic = json.load(infile)
    
    # Initialize empty set of methyls
    methyls = set()
    
    # Map from methyl labels to methyl objects
    label2methyl = {}  

    # Add methyls from the dictionary
    for v in dic["vertices"]:
        
        # Create new methyl from this entry
        new_methyl = Methyl(v["color"], v["seqid"], v["order"],
                            added=v["added"])
        
        # Add to set of new methyls
        methyls.add(new_methyl)
        
        # Create map from label of this methyl to the methyl itself
        label2methyl[new_methyl.label] = new_methyl

    # Create a graph
    structure = nx.Graph()
    structure.add_nodes_from(methyls)

    # Iterate over edges of the graph
    for i, j, distances in dic["edges"]:
        structure.add_edge(label2methyl[i], label2methyl[j],
                           distances=distances)

    # Get set of colors
    colors = [m.color for m in structure.nodes()]
    
    all_colors = set(colors)
    for c in sorted(all_colors):
        print(f"{colors.count(c):2} methyls of type {c}")

    # return the graph
    print()
    return structure


class Methyl:
    """
    Object oriented representation of methyls
    """

    def __init__(self, color, seqid, order, added=False):
        """
        Construct a methyl object of given color, sequence id, and order
        """

        # Save all the features of this residue
        self.color = color  
        self.seqid = seqid
        self.order = order if color in {"L", "V"} else None
        self.label = (f"{color}{seqid}.{order}"
                      if self.order else f"{color}{seqid}")
        self.added = added

    def geminal(self, other):
        """
        Determine if this methyl and another form a geminal pair
        """

        return self.seqid == other.seqid and self.order != other.order

    def to_dict(self):
        """
        Return dictionary representation of self
        """
        return {"color": self.color, "seqid": self.seqid, "order": self.order,
                "label": self.label, "added": self.added}

    def __repr__(self):
        """
        How to represent this as a string
        """

        return f"methyl:{self.label}"

    def __eq__(self, other):
        """
        How to determine if this methyl and another are equal
        """

        return self.label == other.label

    def __hash__(self):
        """
        How to hash this object
        """

        return hash(self.label)


def get_residues(filename: str, model: int, chain: str):
    """
    Read in residues from given filename and model and chain
    """

    # Use BioPython to read in the structure
    models = PARSER.get_structure(f"{filename}", filename)

    # Get the desired model if possible
    if model not in [m.id for m in models]:
        print(f"desired model '{model}' not available in file {filename}")
        print(f"available models: {sorted(m.id for m in models)}")
        exit(1)

    chains = models[model]  # Get all chains in the model

    # Get the desired chain if possible
    if chain not in [c.id for c in chains]:
        print(f"desired chain '{chain}' not available in file {filename}")
        print(f"available chain: {sorted(c.id for c in chains)}")
        exit(1)

    # get the resideues in chain
    return chains[chain].get_residues()  


def get_methyls(filename: str, colors: list, model: int, chain: str) -> list:
    """
    Given the name of a pdb file, return a set of methyl objects extracted
    from the structure
    """

    # Get residues from the given structure model and chain
    residues = get_residues(filename, model, chain)

    # Initialize an empty set of methyls
    methyls = set()  

    # Iterate over residues
    for r in residues:
        
        # get color of this residue
        color = r.get_resname()  
        
        # Ignore this residue if color not desired
        if color not in colors:  
            continue

        # Create two methyls for these
        if color in {"LEU", "VAL"}:  
            methyls.add(Methyl(color[0], r.get_id()[1], 1))
            methyls.add(Methyl(color[0], r.get_id()[1], 2))

        else:
            methyls.add(Methyl(color[0], r.get_id()[1], None))
    
    # return this set of methyls
    return methyls  


def get_atoms(filename: str, colors: list, model: int, chain: str) -> list:
    """
    Return a mapping from methyl labels to atoms
    """

    # Get residues from the given structure model and chain
    residues = get_residues(filename, model, chain)

    # Initialize empty dictionary
    atom_map = {}  

    # Iterate over the residues
    for res in residues:
        
        # get color of this methyl
        color = res.get_resname()  

        # If this res does not have desired color, ignore
        if color not in colors:
            continue

        # Get the sequence id of this residue
        seqid = res.get_id()[1]  
       
        # Get first letter of color
        c = color[0]  

        # Extract the appropriate atoms for methyls of any color
        if color == "LEU":
            atom_map[f"{c}{seqid}.1"] = [a for a in res if 'HD1' in a.id]
            atom_map[f"{c}{seqid}.2"] = [a for a in res if 'HD2' in a.id]

        elif color == "VAL":
            atom_map[f"{c}{seqid}.1"] = [a for a in res if 'HG1' in a.id]
            atom_map[f"{c}{seqid}.2"] = [a for a in res if 'HG2' in a.id]

        elif color == "ALA":
            atom_map[f"{c}{seqid}"] = [a for a in res if 'HB' in a.id]

        elif color == "ILE":
            atom_map[f"{c}{seqid}"] = [a for a in res if 'HD' in a.id]

        elif color == "MET":
            atom_map[f"{c}{seqid}"] = [a for a in res if 'HE' in a.id]

    # return mapping from methyl names to their atoms
    return atom_map  


def pairwise_distance(triplet1, triplet2):
    """
    Takes two triplets of hydrogen atoms and returns the average pairwise
    distances between them
    """

    # Initialize a zero sum
    summation = 0.

    # Iterate over all pairs of hydrogen atoms from this and the other
    for alpha, beta in itertools.product(triplet1, triplet2):

        # Add distance between this pair to the negative sixth
        summation += pow(alpha - beta, -6)

    # Return the average of the sum to the negative 1/6
    return pow(summation/9, -1/6)

