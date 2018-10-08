"""
params.py

Shared parameters for the camera module and the scripts that use it to perform
the methyl assignment problem
"""

# Radius parameters that determine acceptable NOE distance

RADIUS = 10.
ADDED_RADIUS = 10.
SHORT_RADIUS = 10.

# Symmetrization and clustering tolerances

SYM_CTOL = 0.15
SYM_HTOL = 0.02
CLS_CTOL = 0.15
CLS_HTOL = 0.02

# Maximum connected component of symmetrization graph size to use in clustering
# CSP

MAX_COMP_SIZE = 3

# Whether or not assignments or support sets should be forced as hard
# constraints in either CSP

FORCE_ASG = False
FORCE_SV = False

