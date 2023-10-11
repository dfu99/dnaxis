"""
Functions for handling crossover discovery and placement that are specific to SYMMETRIC mode designs
"""

import numpy as np
from ..helper import strandnav
from app import config
from app.routing.motifs import XoverPair


def filtervalid(xovers1, xovers2, threshold):
    valid = []
    for np1 in xovers1:  # match all pairs of potential crossover positions (NuclPairs)
        for np2 in xovers2:
            '''Case: if shape axis is straight'''
            dist = np.linalg.norm(np2.coords - np1.coords)
            if dist < threshold:
                # Just enforce some sanity checks
                if all(n.strand_type() == 'stap' for n in [np1.n5, np1.n3, np2.n5, np2.n3]):
                    # Checks if strands are continuous
                    if all(not strandnav.neargap(n, config.GAP_SPACING_XOVER) for n in [np1.n5, np1.n3, np2.n5, np2.n3]):
                        valid.append(XoverPair(np1, np2))
                    else:  # Routing on strands that are not continuous is not fully supported yet
                        raise NotImplementedError
                elif all(n.strand_type() == 'scaf' for n in [np1.n5, np1.n3, np2.n5, np2.n3]):
                    valid.append(XoverPair(np1, np2))
                else:
                    raise RuntimeError("Nucleotides are on an unknown strand.")
    return valid