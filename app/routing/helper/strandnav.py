# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 21:31:43 2018

@author: Dan
Changelog:
    3/30/2018: Moved helper functions from create.py to this new helper.py file
    4/1/2018: Added angle_between
    4/2/2018: Rename to strandnav
    4/18/2018:
        - Renamed printstruct function to print_to_oxDNA
    7/17/2018:
        - New version of goto5, old_goto5 becomes obsolete
        - New printstrand(base) function for output
    7/27/2018:
        - obsolete functions goto3, goto5, getstrandlength removed
    8/2/2018:
        - added contingency for loops in getstrand
        - added contingency for loops in goto5, goto3
    8/12/2018:
        - Bugfix, getstrand
        - No longer off by 1 due to goto5 moving init_base off by 1
    8/31/2018:
        - getModule modified to accomodate Planes object
    9/6/2018:
        - Changed getModule to getmodule
    11/19/2018:
        added go3, go5 to catch AttributeError when .toThree or .toFive tries to go from a nick
"""


from .. import nucleotide
from . import log


# returns the top level routing that the nucleotide is part of
def getmodule(top, nucl):
    for plane in top:
        for module in plane.modules:
            if nucl in module.structure.scaffold or nucl in module.structure.staple:
                return module


def getregion(base):
    region = []
    # offset so we're not stuck on a 3' feature
    if base.toThree == -1:
        base = base.toFive
    # go to the 5' feature
    # traverse in the 5' direction until reaching a feature
    while not isfeature(base):
        base = base.toFive
    region.append(base)

    base = base.toThree
    region.append(base)
    # traverse in the 3' direction until reaching a feature
    # add all traversed bases to region
    while not isfeature(base):
        base = base.toThree
        region.append(base)
    return region


def checkloop(base):
    start_base = base
    try:
        next_base = base.toThree
    except AttributeError:
        return False  # single base cannot be a loop
    while True:
        if next_base == start_base:
            return True
        elif next_base == -1:
            return False
        else:
            next_base = next_base.toThree


def lenstrand(base):
    """
    Returns length of strand
    :param base: strand, list of Nucleotide
    :return: length, int
    """
    return len(getstrand(base))


def printstrand(base):
    """
    returns the strand sequence
    can start from any base, and will automatically find the 5' end
    20181127 update to use getstrand
    :param base: reference Nucleotide for strand
    :return: string of nucleobases
    """
    sequence = ''
    strand = getstrand(base)
    for n in strand:
        sequence += n.nucl
    return sequence


def getstrand(base):
    """
    Returns the strand as a list of nucleotide objects
    Can start from any base, and will automatically find the 5' end
    :param base: reference Nucleotide for strand
    :return: list of Nucleotide
    """
    sequence = []
    init_base = goto5(base)
    this_base = init_base
    while True:
        sequence.append(this_base)
        try:
            next_base = this_base.toThree
        except AttributeError:
            raise AttributeError("Caution: Strand of length 1.")
            break
        if next_base == -1:
            break
        elif next_base in sequence:  # contingency for loops
            break
        else:
            this_base = next_base
    return sequence


def goto5(init_base):
    """
    from any base, return the base that is at the 5' end of the strand
    input: base is a nucleotide object
    :param init_base: starting Nucleotide(object)
    :return: 5' Nucleotide(object) of strand
    """
    this_base = init_base  # most recently visited
    visited = []
    while True:
        visited.append(this_base)
        next_base = to5(this_base)
        if next_base in visited:  # loop or len 1 contingency
            return init_base
        elif next_base == -1:  # can't go any further
            return this_base
        else:  # continue
            this_base = next_base


def goto3(init_base):
    """
    # does the same thing as goto5 for the other direction
    :param init_base: starting Nucleotide(object)
    :return: 3' Nucleotide(object) of strand
    """
    this_base = init_base  # most recently visited
    visited = []
    while True:
        visited.append(this_base)
        next_base = to3(this_base)
        if next_base in visited:  # loop or len 1 contingency
            return init_base
        elif next_base == -1:  # can't go any further
            return this_base
        else: # continue
            this_base = next_base


def to5(base):
    """
    Moves to 5' Nucleotide with error handler for breaks
    :param base: starting Nucleotide
    :return: adjacent 5' Nucleotide
    """
    try:
        nextbase = base.toFive
    except AttributeError:
        log.system("Can't go any further.")
        return base
    return nextbase


def to3(base):
    """
    Moves to 3' Nucleotide with error handler for breaks
    :param base: starting Nucleotide
    :return: adjacent 3' Nucleotide
    """
    try:
        nextbase = base.toThree
    except AttributeError:
        log.system("Can't go any further.")
        return base
    return nextbase


def absgetstrand(base):
    """
    Returns the non-feature strand as a list of nucleotide objects
    Can start from any base, and will automatically find the 5' end
    :param base: reference Nucleotide for strand
    :return: list of Nucleotide
    """
    sequence = []
    init_base = absgoto5(base)
    this_base = init_base
    while True:
        sequence.append(this_base)
        try:
            next_base = this_base.__strand3__
        except AttributeError:
            log.system("Caution: Strand of length 1.")
            raise AttributeError
            break
        if next_base == -1:
            break
        elif next_base in sequence:  # contingency for loops
            break
        else:
            this_base = next_base
    return sequence


def absgoto5(init_base):
    """
    from any base, return the base that is at the 5' end of the strand using absolute adjacency
    input: base is a nucleotide object
    :param init_base: starting Nucleotide(object)
    :return: 5' Nucleotide(object) of strand
    """
    this_base = init_base  # most recently visited
    visited = []
    while True:
        visited.append(this_base)
        next_base = absto5(this_base)
        if next_base in visited:  # loop or len 1 contingency
            return init_base
        elif next_base == -1:  # can't go any further
            return this_base
        else:  # continue
            this_base = next_base


def absgoto3(init_base):
    """
    # does the same thing as absgoto5 for the other direction
    :param init_base: starting Nucleotide(object)
    :return: 3' Nucleotide(object) of strand
    """
    this_base = init_base  # most recently visited
    visited = []
    while True:
        visited.append(this_base)
        next_base = absto3(this_base)
        if next_base in visited:  # loop or len 1 contingency
            return init_base
        elif next_base == -1:  # can't go any further
            return this_base
        else:  # continue
            this_base = next_base


def absto5(base):
    """
    Moves to 5' Nucleotide with error handler for breaks using absolute adjacency
    :param base: starting Nucleotide
    :return: adjacent 5' Nucleotide
    """
    try:
        nextbase = base.__strand5__
    except AttributeError:
        log.system("Can't go any further.")
        return base
    return nextbase


def absto3(base):
    """
    Moves to 3' Nucleotide with error handler for breaks using absolute adjacency
    :param base: starting Nucleotide
    :return: adjacent 3' Nucleotide
    """
    try:
        nextbase = base.__strand3__
    except AttributeError:
        log.system("Can't go any further.")
        return base
    return nextbase


def get5strand(strand):
    """
    Returns the 5' side strand
    Input 'strand' is 5'-3'
    :param strand: list of Nucleotides of current strand
    :return: list of Nucleotides of 5' adjacent strand
    """
    n = strand[0]
    n5 = n.__strand5__
    adj5 = getstrand(n5)
    return adj5


def get3strand(strand):
    """
    Returns the 3' side strand
    Input 'strand' is 5'-3'
    :param strand: list of Nucleotides of current strand
    :return: list of Nucleotides of 3' adjacent strand
    """
    n = strand[-1]
    n3 = n.__strand3__
    adj3 = getstrand(n3)
    return adj3


def compare_min(s1, s2):
    """
    returns the shortest of 2 strands
    :param s1: strand 1, list of Nucleotide
    :param s2: strand 2, list of Nucleotide
    :return: strand, list of Nucleotide
    """
    if len(s1) < len(s2):
        return s1
    return s2


def strandequals(strand1, strand2):
    """
    Compares the sequences of two strands with permutations accepted.
    :param strand1: list of Nucleotide
    :param strand2: list of Nucleotide
    :return: True/False = Equal/Unequal
    """
    if all(s1 in strand2 for s1 in strand1) and all (s2 in strand1 for s2 in strand2):
        return True
    return False


def seqequals(strand1, strand2):
    """
    Compares the sequences of two strands with no permutations.
    :param strand1: list of Nucleotide
    :param strand2: list of Nucleotide
    :return: True/False = Equal/Unequal
    """
    '''Error handlers'''
    if type(strand1) != list or type(strand2) != list:
        # Not comparing lists
        return False
    if any(type(n) != nucleotide.Nucleotide for n in strand1) or any(type(n) != nucleotide.Nucleotide for n in strand2):
        # Not comparing lists of nucleotides
        return False
    ''''''
    for s1, s2 in zip(strand1, strand2):
        if s1 != s2:
            return False
    return True


def distfromfeature(nucl):
    """
    Calculates the bp distance to the nearest feature (nick or crossover)
    :param nucl: Nucleotide object denoting the reference position
    :return: int distance in bp counted to next feature
    """
    d = 0
    if isfeature(nucl) or isfeature(nucl.Comp):
        return d
    n5 = nucl.__strand5__
    n3 = nucl.__strand3__
    n5c = n5.Comp
    n3c = n3.Comp
    while not any((isfeature(n5), isfeature(n3), isfeature(n5c), isfeature(n3c), n5 == nucl, n3 == nucl)):
        n5 = n5.__strand5__
        n3 = n3.__strand3__
        n5c = n5c.toFive
        n3c = n3c.toThree
        d += 1
    return d


def distfromxover(nucl, options="both"):
    """
    Looks for distance to crossovers nearest to the input reference position
    :param nucl: Nucleotide object reference position
    :return: int distance in bp counted to next crossover
    """
    d = 0
    n5 = nucl.__strand5__
    n3 = nucl.__strand3__
    n5c = n5.Comp
    n3c = n3.Comp
    if isxover(nucl) or isxover(nucl.Comp):
        return d
    elif options == "this":
        check = (isxover(n5), isxover(n3))
    elif options == "comp":
        check = (isxover(n5c), isxover(n3c))
    elif options == "both":
        check = (isxover(n5), isxover(n3), isxover(n5c), isxover(n3c))
    else:
        raise RuntimeError("distfromxover does not support option {}".format(options))

    while not any(check):
        n5 = n5.__strand5__
        n3 = n3.__strand3__
        n5c = n5c.__strand5__
        n3c = n3c.__strand3__
        d += 1
        if options == "this":
            check = (isxover(n5), isxover(n3))
        elif options == "comp":
            check = (isxover(n5c), isxover(n3c))
        elif options == "both":
            check = (isxover(n5), isxover(n3), isxover(n5c), isxover(n3c))
        else:
            raise RuntimeError("distfromxover does not support option {}".format(options))
    return d


def xp_distfromxover(xp):
    """
    Looks for distance to crossovers nearest to the input reference position
    :param xp: Initial position at XoverPair object
    :return: int distance in bp counted to next crossover
    """
    d = 1
    init_nucl = [xp.n1.n5, xp.n1.n3, xp.n2.n5, xp.n2.n3]
    check_nucl = splash(*init_nucl)
    check_nuclc = [n.Comp for n in check_nucl]
    while True:
        if any(isxover(n) for n in check_nucl + check_nuclc):
            return d
        elif any(n in init_nucl for n in check_nucl):
            return d
        else:
            check_nucl = splash(*check_nucl)
            check_nuclc = [n.Comp for n in check_nucl]
            d += 1


def xover_common_edge(x1, x2, same_edge):
    """
    Checks if two XoverPairs span the same modules (are the same edge)
    :param x1: XoverPair
    :param x2: XoverPair
    :param same_edge: Must exist on same edge (True) or can also exist on adjacent edges (False)
    :return: True/False
    """
    x1n1 = x1.n1.n5
    x1n2 = x1.n2.n3
    x2n1 = x2.n1.n5
    x2n2 = x2.n2.n3

    x1m1 = x1n1.get_top_module()
    x1m2 = x1n2.get_top_module()
    x2m1 = x2n1.get_top_module()
    x2m2 = x2n2.get_top_module()
    if same_edge:  # Checks same edge
        if all(m in [x2m1, x2m2] for m in [x1m1, x1m2]) and all(m in [x1m1, x1m2] for m in [x2m1, x2m2]):
            return True
    else:  # Checks any adjacent edge
        if any(m in [x2m1, x2m2] for m in [x1m1, x1m2]) and any(m in [x1m1, x1m2] for m in [x2m1, x2m2]):
            return True
    return False


def splash(n1n5, n1n3, n2n5, n2n3):
    """
    Iterate outwards from a 4-tuple originating from a crossover
    :param n1n5:
    :param n1n3:
    :param n2n5:
    :param n2n3:
    :return:
    """
    return [n1n5.__strand5__, n1n3.__strand3__, n2n5.__strand5__, n2n3.__strand3__]


def neargap(nucl, limit):
    """
    Checks whether nucleotide is far enough from a acyclic gap
    Returns True/False rather than distance because unnecessary to search past certain distance
    :param nucl: Nucleotide object
    :param limit: The search distance
    :return: True/False
    """
    this_nucl = nucl
    if this_nucl.is_on_gap():
        return True
    n5 = n3 = this_nucl
    for i in range(limit):
        n5 = n5.__strand5__
        n3 = n3.__strand3__
        if n5.is_on_gap() or n3.is_on_gap():
            return True
    return False


def isfeature(nucl):
    """
    Checks whether a Nucleotide is a nick or crossover
    :param nucl: Nucleotide object to check
    :return: True/False
    """
    if isnick(nucl):
        return True
    elif isxover(nucl):
        return True
    return False


def isxover(nucl):
    """
    Checks whether a Nucleotide is a crossover
    :param nucl: Nucleotide object to check
    :return: True/False
    """
    try:
        nucl_helix = nucl.get_top_strand()
        if not (nucl.toThree.get_top_strand() == nucl_helix and
                nucl.toFive.get_top_strand() == nucl_helix):
            return True
        return False
    except AttributeError:
        # If nucl is a nick then toThree or toFive will == -1
        return False


def nucl_same_helix(n1, n2):
    """
    Checks if both nucls can be reached in the same helix
    :param n1: Nucleotide object
    :param n2: Nucleotide object
    :return: True/False
    """
    strand = absgetstrand(n1)
    if n2 in strand:
        return True
    else:
        return False


def isXover5End(nucl):
    if isxover(nucl) and isxover(nucl.toThree):
        return True
    return False


def isXover3End(nucl):
    if isxover(nucl) and isxover(nucl.toFive):
        return True
    return False


def isnick(nucl):
    """
    Checks whether a Nucleotide is a nick
    :param nucl: Nucleotide object to check
    :return: True/False
    """
    if nucl.toThree == -1 or nucl.toFive == -1:
        return True
    return False
