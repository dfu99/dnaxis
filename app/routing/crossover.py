#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  8 20:21:24 2018

@author: dfu

crossover.py

CHANGELOG:
    9/8/2018:
        File created. Separated crossover functions from origami.py
    9/12/2018:
        is_xover3, is_xover5, and get_xover_complement added
        Added top level functions from sym_origami.py
    9/29/2018
        Create 2 new classes, FullXover and HalfXover
    10/08/2018
        Added onmodule to FullXover
    10/10/2018
        Added XoverPair class
        Added NuclPair class
    10/11/2018
        Removed obsolete functions
        Forked strand specific methods to stxover.py, scxover.py
        
"""
from .helper import log, strandnav
from .mode_symmetric import crossover as symxovers
from .mode_asymmetric import crossover as asymxovers

from . import solve_crossovers

from app import config
from .exceptions import *

from operator import itemgetter

import itertools as itls
import math

class Mixin:
    # finds number of bases traversed until reaching an xover
    # TODO: cnt = 1 is not accurate for cases starting on an xover
    def dist2nearestxover(self, nucl):
        nucl3 = nucl.toThree
        nucl5 = nucl.toFive
        cnt = 1
        while True:
            if self.is_xover(nucl3) or self.is_xover(nucl5):
                return cnt
            else:
                nucl3 = nucl3.toThree
                nucl5 = nucl5.toFive
                cnt += 1

    def clear_seam(self):
        """
        Looks for a seam ands removes the violating staple crossover
        :return: None
        """
        all_objxovers = self.get_all_xovers()
        for xover in all_objxovers:
            if xover.strandtype == 'scaf':
                (n1, n2) = (xover.n1, xover.n2)
                (n1_5, n1_3) = n1.to_tuple
                (n2_5, n2_3) = n2.to_tuple

                for i in range(4):
                    n1_5 = n1_5.__strand5__
                    n1_3 = n1_3.__strand3__
                    n2_5 = n2_5.__strand5__
                    n2_3 = n2_3.__strand3__
                    n1_5c = n1_5.Comp
                    n1_3c = n1_3.Comp
                    n2_5c = n2_5.Comp
                    n2_3c = n2_3.Comp
                    check_list = [n1_5c, n1_3c, n2_5c, n2_3c]
                    for c in check_list:
                        if strandnav.isxover(c):
                            unxover = self.findxover(c)
                            log.system("A seam was found and the staple crossover was removed between modules "
                                       "{} and {}: {}".format(n1.up(), n2.up(), xover))
                            unxover.undo()

    def apply_crossovers(self, xoverset, strandtype):
        """
        apply all crossovers from the xoverset
        :param xoverset: list of motifs.XoverPair
        :param strandtype: 'scaf' or 'stap'
        :return: None
        """
        if strandtype == 'scaf':
            protect_val = config.SCPROTECT
        elif strandtype == 'stap':
            protect_val = config.PROTECT
        else:
            raise TypeError("apply_crossovers() got an unexpected stand type '{}'".format(strandtype))
        for xp in xoverset:
            fx = xp.apply()
            self.protect(xp.n1.to_tuple[0], protect_val)
            self.protect(xp.n2.to_tuple[0], protect_val)
            fx.settype(strandtype)
            xp.n1.up().addxover(fx)
            xp.n2.up().addxover(fx)

# =============================================================================
# Static functions
# =============================================================================


def get_gap_xovers(module1, module2):
    """
    Creates a special gap spanning crossover for pairs of adjacent non-cyclic modules
    :param module1:
    :param module2:
    :return:
    """

    # print("DEBUG: Xovers to choose from: ", xoverset)
    # Find the gap
    # Save the gap endpoints as gap5 and gap3
    for nucl in module1.helix.scaf:
        if nucl.__strand5__ == -1:
            m1gap3 = nucl
        if nucl.__strand3__ == -1:
            m1gap5 = nucl

    for nucl in module2.helix.scaf:
        if nucl.__strand5__ == -1:
            m2gap3 = nucl
        if nucl.__strand3__ == -1:
            m2gap5 = nucl

    # From the two half crossovers, build another set of NuclPairs
    # @invariant, 1st index is the nucleotide (this_nucl) used during traversal
    # Traversing from 5', the nucleotide is definitely the 3' of the NuclPair
    n1_3 = m1gap3
    #   Then the other half of the HalfXover is the 5' of the opposing NuclPair
    n2_5 = m2gap5
    # Traversing from 3', the nucleotide is definitely the 5' of the NuclPair
    n1_5 = m1gap5
    #   Then the other half of the HalfXover is the 3' of the opposing NuclPair
    n2_3 = m2gap3
    npair1 = GapNuclPair(n1_5, n1_3)
    # TODO: Get rid of this, this is HORRIBLE
    n1_5.get_top_origami().gapnuclpairs.append(npair1)
    # TODO: Get rid of this, this is HORRIBLE
    npair2 = GapNuclPair(n2_5, n2_3)
    n1_5.get_top_origami().gapnuclpairs.append(npair2)
    # Return the XoverPair
    return [XoverPair(npair1, npair2)]


def filtertoregion(region, xoverpairs):
    """
    Keep only XoverPairs that exist within the specified strand
    :param region: list of Nucleotides
    :param xoverpairs: list of XoverPairs
    :return: reduced list of XoverPairs
    """
    valid = []
    for xp in xoverpairs:
        if any(n in xp for n in region):
            valid.append(xp)
    return valid


def filtervalid(xovers1, xovers2, routing_mode="symmetric", **kwargs):
    """
    From a pair of lists of NuclPairs indicating valid crossover positions
    find any pair that has acceptable separation (< VALIDXOVERTHRESHBP)

    :param xovers1: list of NuclPair
    :param xovers2: list of NuclPair
    :param routing_mode: asymmetric or symmetric
    :return: List of XoverPair objects
    """
    threshold = kwargs.get("thresh", config.VALIDXOVERTHRESHBP)
    valid = []

    # Quick error check
    if not xovers1 or not xovers2:
        # raise NoValidXoverThreshold("Potential alignment positions are empty.")
        log.log_warning("Potential alignment positions are empty.")
        return valid

    if routing_mode == "asymmetric":
        valid = asymxovers.filtervalid(xovers1, xovers2, threshold)
    elif routing_mode == "symmetric":
        valid = symxovers.filtervalid(xovers1, xovers2, threshold)
    else:
        raise NotImplementedError
    # log.log_debug("Returning [filtervalid] with {} crossovers".format(len(valid)))
    return valid


def filteroverlap(xoverpairs):
    """
    Check for any xoverpairs that are already created crossovers to avoid overlap.
    :param xoverpairs: distance sorted list of XoverPair
    :return: reduced list of XoverPair
    """
    valid = []
    for xp in xoverpairs:
        (xp1n5, xp1n3) = xp.n1.to_tuple
        (xp2n5, xp2n3) = xp.n2.to_tuple
        if any(strandnav.isfeature(nucl) for nucl in [xp1n5, xp1n3, xp2n5, xp2n3]):
            pass
        else:
            valid.append(xp)
    return valid


def filterdist2feat(xoverpairs,
                    spacing_init=config.VALIDXOVERSPACING_SAME,
                    spacing_limit=config.LIMSPACINGOFFSET,
                    suppress_novalidspacing=False, **kwargs):
    """
    Iterates through a list of XoverPair objects.
    Filters for and then sorts by pairs that are farthest from any existing crossover.
    Only detects FullXovers
    :param xoverpairs: list of XoverPairs
    :param spacing_init: initial space constraint
    :param spacing_limit: iterative cap
    :param suppress_novalidspacing: suppress errors
    :return: reduced list of XoverPairs
    """
    valid = []
    spacingoffset = 0
    numxovers = kwargs.get("numxovers", 0)
    while True:
        for xp in xoverpairs:
            mindist = strandnav.xp_distfromxover(xp)
            if mindist > spacing_init - spacingoffset:
                valid.append((xp, mindist))
        if (not valid) or len(valid) < numxovers:
            spacingoffset += 1
            if spacingoffset > spacing_limit and not suppress_novalidspacing:
                raise NoValidSpacing("[filterdist2feat]"
                                     "No satisfying space for additional crossovers was found. Consider decreasing "
                                     "the initial base spacing ({}) or increasing the spacing relax limit ({})."
                                     .format(spacing_init, spacing_limit))
            else:
                continue
        else:
            break
    sortvalid = sorted(valid, key=itemgetter(1), reverse=True)
    res = [tpl[0] for tpl in sortvalid]
    return res  # return the xoverpairs in distance sorted order


def distribute(xoverpairs, numxovers, strandtype):
    """
    Distributes a fixed number of crossovers as evenly as possible
    :param xoverpairs: list of XoverPair
    :param numxovers: int value
    :return:reduced set of XoverPair
    """
    arr, ref_module = xovers_to_arr(xoverpairs, strandtype)
    if numxovers == 1:
        s_arr = ([arr[math.ceil(len(arr)/2)]], 0.0, True)
    else:
        s_arr = solve_crossovers.largestMinDist(arr, numxovers)

    # Flag error if solve_crossovers didn't find a solution
    if not s_arr[2]:
        raise RuntimeError("No valid distribution of {} XoverPairs was found.".format(numxovers))

    # Consolidate the reference strand
    if strandtype == 'scaf':
        ref_strand = ref_module.helix.scaf
    elif strandtype == 'stap':
        ref_strand = ref_module.helix.stap
    else:
        raise KeyError("Incompatible strand type argument '{}'".format(strandtype))

    # Convert back to XoverPairs
    x_arr = arr_to_xovers(s_arr[0], xoverpairs, ref_strand)

    # Return the chosen XoverPairs
    return x_arr


def filterspacing(xoverpairs):
    """
    Selects a set of XoverPairs from the provided list that fulfill the configured spacing between each adjacent
    crossover.
    :param xoverpairs: List of XoverPair (sorted sequentially by position)
    :return: Reduced list of XoverPair
    """
    # Convert the crossovers to integer positions
    arr, ref_module = xovers_to_arr(xoverpairs, 'stap')

    # Solve for the set of values satisfying the spacing
    if len(arr) == 1:
        s_arr = (arr, 0.0)
    else:
        # print("Looking for largest min {} distance in {}".format(minspace, arr))
        s_arr = solve_crossovers.iterate(arr, config.XMDYN_MINSPACEOFF)

    # Return no crossovers if solve_crossovers didn't find a solution
    if s_arr[1] == -1:
        return []

    # Dynamic xover mode only matters for staple crossovers
    ref_strand = ref_module.helix.stap

    # Convert back to XoverPairs
    x_arr = arr_to_xovers(s_arr[0], xoverpairs, ref_strand)

    # Return the chosen XoverPairs
    return x_arr


def filterduplicate(xoverpairs):
    """
    Check for any XoverPairs that may use the same NuclPair and thus collide when they are applied
    :param xoverpairs: distance sorted list of XoverPairs
    :return: reduced list of XoverPairs
    """
    for a, b in itls.combinations(xoverpairs, 2):
        xp1 = a
        xp2 = b

        (xp1_n1_5, xp1_n1_3) = xp1.n1.to_tuple
        (xp1_n2_5, xp1_n2_3) = xp1.n2.to_tuple
        xp1_nucls = (xp1_n1_5, xp1_n1_3, xp1_n2_5, xp1_n2_3)
        (xp2_n1_5, xp2_n1_3) = xp2.n1.to_tuple
        (xp2_n2_5, xp2_n2_3) = xp2.n2.to_tuple
        xp2_nucls = (xp2_n1_5, xp2_n1_3, xp2_n2_5, xp2_n2_3)

        for n1 in xp1_nucls:
            for n2 in xp2_nucls:
                if n1 is n2:
                    xp1_str = [n.numid for n in xp1_nucls]
                    xp2_str = [n.numid for n in xp2_nucls]
                    log.debug("{} uses the same nucleotides as {}".format(xp1_str, xp2_str))
                    try:
                        xoverpairs.remove(b)
                    except ValueError:
                        pass
    return xoverpairs


def sort_xovers(xoverpairs):
    """
    Sort the list of XoverPair object in sequential order
    :param xoverpairs: List of XoverPair
    :return: sorted list of XoverPair
    """
    # Use the tighter module as reference
    ref_xover = xoverpairs[0]
    m1 = ref_xover.n1.up()
    m2 = ref_xover.n2.up()
    if m1.bp <= m2.bp:
        ref_module = m1
    else:
        ref_module = m2
    # Dynamic xover mode only matters for staple crossovers
    ref_strand = ref_module.helix.stap
    # Initialize output
    ordered_xp = []

    for this_nucl in ref_strand:
        # Check all the XoverPairs for this nucleotide
        for xp in xoverpairs:
            # Check by their 5' endpoint only to avoid double counting
            if xp.n1.n5 == this_nucl or xp.n2.n5 == this_nucl:
                ordered_xp.append(xp)

    # If loop exits, we've gone around the entire ring
    # Check that all the XoverPairs were accounted for
    for xp in xoverpairs:
        if xp in ordered_xp:
            pass
        else:
            raise RuntimeError("Sorted crossover pairs do not match the input.")

    return ordered_xp


# nextxovers are potential xovers to a desired ring
# find the xover in nextxovers that is the greatest distance from any existing xover
def get_bestspace(nextxovers):
    d = 0
    for nx in nextxovers:
        (n5, n3) = nx.to_tuple
        _d = min(strandnav.distfromfeature(n5), strandnav.distfromfeature(n3))
        if _d > d:
            d = _d
            bestx = nx
    if d < config.VALIDXOVERSPACING_SAME:
        thisring = nextxovers[0].up()
        log.out(__name__, "Warning! There's no more space on {} to place crossovers.".
                format(thisring.getsimpleposition()) +
                "The best option is {} with {} distance to closest feature.".
                format(nx.numid, d))
    return bestx


def filter_anti_parallel(xoverpairs):
    """
    Remove parallel crossovers
    :param xoverpairs: List of XoverPair potential crossover positions
    :return: List
    """
    filtered_list = []
    for xp in xoverpairs:
        if xp.is_anti_parallel():
            filtered_list.append(xp)
    return filtered_list


def get_xovers(module1, module2, strand, **kwargs):
    """
    Get the number (num) of valid (filtervalid) and correctly spaced (filterdist2feat) crossovers
    :param module1: source connecting module
    :param module2: target connecting module
    :param strand: options: staple ('stap')
                            scaffold ('scaf')
    :return: list of XoverPair objects to later apply
    """

    xovers1 = module1.getxoversto(module2, strand)
    xovers2 = module2.getxoversto(module1, strand)
    if not xovers1 or not xovers2:
        raise RuntimeError("No crossovers were found between modules {} and {}. Please confirm they are adjacent.\n"
                           "xovers1: {}\n"
                           "xovers2: {}\n".format(module1, module2, xovers1, xovers2))
    bondlen = kwargs.get("threshold", config.VALIDXOVERTHRESHBP)
    num = kwargs.get('numxovers', 0)
    while True:
        try:
            # Filter by bond length
            xplist = filtervalid(xovers1, xovers2, thresh=bondlen)
            log.debug(__name__, "Num Valid Xovers: {}".format(len(xplist)))
            # Filter positions overlapping with existing xovers
            xplist = filteroverlap(xplist)
            log.debug(__name__, "Num Non-overlapping Xovers: {}".format(len(xplist)))
            # Filter positions too close to existing xovers
            xplist = filterdist2feat(xplist, numxovers=num)
            log.debug(__name__, "Num Distance Xovers: {}".format(len(xplist)))
            # Filter duplicates in the list
            xplist = filterduplicate(xplist)
            log.debug(__name__, "Num Duplicate Xovers: {}".format(len(xplist)))
            if strand == 'scaf':  # Enforce anti-parallel scaffold crossovers
                xplist = filter_anti_parallel(xplist)
                log.debug(__name__, "Num Anti-Parallel Xovers: {}".format(len(xplist)))

        # Suppress errors for empty sets if XOVERMODE dynamic
        except NoValidSpacing:
            if config.XOVERMODE == 'dynamic' and strand == 'stap':
                return []
            else:  # Pass through error
                log.system("WARNING: Increasing bond length locally from {:2.1f} to {:2.1f}".format(bondlen, bondlen+0.1))
                bondlen = round(bondlen + 0.1, 1)
                continue
        except NoValidXoverThreshold:
            if config.XOVERMODE == 'dynamic' and strand == 'stap':  # DEPRECATED?
                return []
            else:  # Pass through error
                log.system("WARNING: Increasing bond length locally from {:2.1f} to {:2.1f}".format(bondlen, bondlen+0.1))
                bondlen = round(bondlen + 0.1, 1)
                continue

        # Error flagging
        if not xplist:  # Empty
            log.log_warning("The crossover set between \n\n\t{}#{}\n\n\tand\n\n\t{}#{}\n\tis empty."
                            .format(str(module1), module1.up().up().get_module_index(module1),
                                    str(module2), module2.up().up().get_module_index(module2)))
        if 'numxovers' in kwargs:
            if len(xplist) < num:
                log.log_warning("The crossover set (size {}) between "
                                "\n\n\t{}\n\n\tand\n\n\t{}\n\t"
                                "does not satisfy the requested number ({}) of crossovers."
                                .format(len(xplist), str(module1), str(module2), num))

        return xplist


def xovers_to_arr(xoverpairs, strandtype):
    """
    Converts a list of XoverPair potential xover positions to a integer array
    :param xovers: Sorted List of XoverPair
    :return: Ordered list of integer positions
    """

    # Use the tighter module as reference
    ref_xover = xoverpairs[0]
    m1 = ref_xover.n1.up()
    m2 = ref_xover.n2.up()
    if m1.bp <= m2.bp:
        ref_module = m1
    else:
        ref_module = m2
    if strandtype == 'stap':
        ref_strand = ref_module.helix.stap
    elif strandtype == 'scaf':
        ref_strand = ref_module.helix.scaf
    else:
        raise KeyError("Incompatible strand type argument '{}'".format(strandtype))

    # Initialize output
    poslist = []

    for pos, this_nucl in enumerate(ref_strand, start=1):
        # Check all the XoverPairs for this nucleotide
        for xp in xoverpairs:
            # Check by their 5' endpoint only to avoid double counting
            if xp.n1.n5 == this_nucl or xp.n2.n5 == this_nucl:
                poslist.append(pos)

    # Check that all the XoverPairs were accounted for
    if len(poslist) != len(xoverpairs):
        raise RuntimeError("Position list of size {} does not match input of size {}.".
                           format(len(poslist), len(xoverpairs)))

    if not poslist:
        raise RuntimeError("XoverPair to array conversion is empty.")

    return poslist, ref_module


def arr_to_xovers(arr, xoverpairs, ref_strand):
    """
    Backwards conversion from array indices back to crossovers from the given xoverpairs list
    :param arr: array
    :param xoverpairs: list of possible XoverPairs
    :param ref_strand: list of nucleotides (strand) of helix
    :return:
    """
    # Look for each position
    iter_arr = iter(arr)
    findpos = next(iter_arr)
    # print("Looking for", findpos)
    # Convert those values back to the corresponding XoverPair objects
    x_arr = []
    for pos, this_nucl in enumerate(ref_strand, start=1):
        try:
            # Check all the XoverPairs for this nucleotide
            if pos == findpos:
                for xp in xoverpairs:
                    # Check by their 5' endpoint only to avoid double counting
                    if xp.n1.n5 == this_nucl or xp.n2.n5 == this_nucl:
                        x_arr.append(xp)
                findpos = next(iter_arr)
        except StopIteration:
            break
    return x_arr


def sortdist2feat(xoverpairs, param="both"):
    """
    Iterates through a list of XoverPair objects.
    Sorts by pairs that are farthest from any existing crossover.
    :param xoverpairs:
    :return:
    """
    valid = []
    for xp in xoverpairs:
        (xp1n5, xp1n3) = xp.n1.to_tuple
        (xp2n5, xp2n3) = xp.n2.to_tuple
        points = (xp1n5, xp1n3, xp2n5, xp2n3)
        if param == "both":
            mindist = min([strandnav.distfromfeature(n) for n in points])
        elif param == "xover":
            mindist = min([strandnav.distfromxover(n) for n in points])
        valid.append((xp, mindist))

    sortvalid = sorted(valid, key=itemgetter(1), reverse=True)
    if sortvalid[0][1] < 0:
        # print("Even the best choice is too close to an existing crossover. That's going to be a problem.")
        return []
    res = [tpl[0] for tpl in sortvalid]
    return res  # return the xoverpairs in distance sorted order
