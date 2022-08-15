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
from .helper import dnaconnector, mymath, mytrig, log, rmatrix, strandnav

from . import nucleotide, modules, solve_crossovers

from app import config
from .exceptions import *

from operator import itemgetter
import numpy as np
import itertools as itls
import math


class Xover:
    """
    General Xover methods
    """
    def is_anti_parallel(self):
        """
        Checks if this crossover is a parallel or anti-parallel crossover
        :return: True/False
        """
        # Only need to check one half crossover
        # Dot product the direction vectors.
        # If the difference is more than 180 degrees then the strands are anti parallel
        # If the two half xovers are at a large separation distance, then check if the strands are also acyclic
        #   and if so, ignore because this is a gap-spanning crossover that will probably return a false negative
        m1 = self.n1.n5.get_top_module()
        m2 = self.n2.n5.get_top_module()
        if np.linalg.norm(self.n1.n5.center - self.n1.n3.center) > 3:  # and noncyclic.gap_exists(m1, m2):
            return True
        d1 = self.n1.n5.duvec
        d2 = self.n2.n3.duvec
        dp = np.dot(d1, d2)
        t = np.degrees(np.arccos(dp / (np.linalg.norm(d1) * np.linalg.norm(d2))))
        if t >= 90:
            return True
        else:
            # highlight_nucl(self.n1.n5)
            return False


class FullXover(Xover):
    """
    Forms a full crossover connection between pairs of nucleotide pairs
    Input two nucleotides (each is the 5' base) of two adjacent strands
    Assumes n1_5 and n2_5 are diagonal positions on adjacent anti-parallel helices
    ======53=======
    ======35=======
    ======53=======
    A full crossover is 2 half crossovers
    """
    def __init__(self, n1, n2):
        self.n1 = n1
        self.n2 = n2
        (n1_5, n1_3) = n1.to_tuple
        (n2_5, n2_3) = n2.to_tuple
        self.xovers = (HalfXover(n1_5, n2_3), HalfXover(n2_5, n1_3))
    
    # restores old connections
    def undo(self):
        (n1_5, n1_3) = self.n1.to_tuple
        (n2_5, n2_3) = self.n2.to_tuple
        # get the rings
        r1 = self.n1.up()
        r2 = self.n2.up()
        # remove the crossover from the rings
        r1.delxover(self)
        r2.delxover(self)
        # remove reference to HalfXovers
        self.xovers = None
        # reconnect original connections
        log.debug(__name__, "Undo crossover", self.n1.numid, "-", self.n2.numid)
        dnaconnector.nuclconnect(n1_5, n1_3)
        dnaconnector.nuclconnect(n2_5, n2_3)

    def onmodule(self, module):
        """
        The crossover does not have any built-in directionality, thus we do not know which pair of nucleotides
            is on which side of the crossover (and on which strand)
        Helps discern which NuclPair is on which module
        :param module: the module object
        :return: the NuclPair that is on module
        """
        m1 = self.n1.up()
        m2 = self.n2.up()
        if m1 == module:
            return self.n1
        elif m2 == module:
            return self.n2
        else:
            raise ValueError("This crossover does not exist on Ring({},{}).".format(module.bp, module.height))

    def __contains__(self, value):
        """
        Checks if Nucleotide value is present in the XoverPair
        :param value:
        :return:
        """
        if value in (self.n1.n5, self.n1.n3, self.n2.n5, self.n2.n3):
            return True
        else:
            return False

    def spans(self, m1, m2):
        """
        Returns True/False whether crossover spans the input strands
        :param m1: Any shape or strand object
        :param m2: Any shape or strand object
        :return: True/False
        """
        if issubclass(type(m1), modules.Shape) and issubclass(type(m2), modules.Shape):
            nt1 = self.n1.n5.get_top_module()
            nt2 = self.n2.n5.get_top_module()

        elif type(m1) == modules.JoinedStrand and type(m2) == modules.JoinedStrand:
            nt1 = self.n1.n5.get_top_strand()
            nt2 = self.n2.n5.get_top_strand()

        # Checks equality by their topology value (should be correctly set)
        if all(x in [m1, m2] for x in [nt1, nt2]):
            return True
        else:
            return False

    def has_nucl(self, nucl):
        """
        Returns True/False whether Nucleotide nucl comprises the FullXover or not
        :param nucl: Nucleotide
        :return: True/False
        """
        if nucl in [self.n1.n5, self.n1.n3, self.n2.n5, self.n2.n3]:
            return True
        else:
            return False

    def settype(self, choose_strand):
        if choose_strand not in ['scaf', 'stap']:
            raise ValueError("Strand type must be 'scaf' or 'stap'.")
        self.strandtype = choose_strand

    def __repr__(self):
        h1 = self.xovers[0]
        h2 = self.xovers[1]
        return "\n{:^15} {:^15}\n{:^15} {:^15}\n{:^15} {:^15}".format(
            h1.nucls[0].numid, h2.nucls[1].numid, "|", "|", h1.nucls[1].numid, h2.nucls[0].numid)

    def __eq__(self, other):
        """
        Checks equivalence by matching all 4 nucleotide crossover points
        :param other: FullXover
        :return: True/False
        """
        if type(other) != type(self):
            return False
        x1n1n5 = self.n1.n5
        x1n1n3 = self.n1.n3
        x1n2n5 = self.n2.n5
        x1n2n3 = self.n2.n3
        x1 = [x1n1n5, x1n1n3, x1n2n5, x1n2n3]
        x2n1n5 = other.n1.n5
        x2n1n3 = other.n1.n3
        x2n2n5 = other.n2.n5
        x2n2n3 = other.n2.n3
        x2 = [x2n1n5, x2n1n3, x2n2n5, x2n2n3]
        if all(nid in x2 for nid in x1) and all(nid in x1 for nid in x2):
            return True
        else:
            return False

    def complement(self):
        """
        Returns a virtual complement as an XoverPair
        :return: XoverPair
        """
        (n1_5, n1_3) = self.n1.to_tuple
        (n2_5, n2_3) = self.n2.to_tuple
        n1 = NuclPair(n1_3.Comp, n1_5.Comp, nicked=True)
        n2 = NuclPair(n2_3.Comp, n2_5.Comp, nicked=True)
        xp = XoverPair(n1, n2)
        return xp


class HalfXover:
    def __init__(self, n_5, n_3):
        self.nucls = (n_5, n_3)
        dnaconnector.nuclconnect(n_5, n_3)

    def __repr__(self):
        return "{}-{}".format(self.nucls[0].numid, self.nucls[1].numid)


class XoverPair(Xover):
    """
    facilitate the format of an xover pair as type tuple
    define: an xover pair is a pair of nucleotide pairs on adjacent helices
      where an xover could be placed
    an xover pair should be formatted:
    (<nucleotide pair (5',3') #1>,
    <nucleotide pair (5',3') #2>,
    <distance between the nucleotide pairs>,
    <angular position of nucleotide pair #1>,
    <angular position of nucleotide pair #2>)
    """
    def __init__(self, n1, n2):
        self.n1 = n1  # n1 is NuclPair object
        self.n2 = n2  # n2 is NuclPair object
        self.dist = abs(mytrig.angle_diff(self.n1.pos, self.n2.pos))
        
    def fget(self):
        return self.n1, self.n2, self.dist, self.n1.pos, self.n2.pos
    
    def fset(self, *value):
        (self.n1, self.n2, self.dist) = value
        
    def apply(self):
        return FullXover(self.n1, self.n2)

    def get_halfxover(self, nucl):
        """
        Gets the HalfXover nucleotides corresponding to the input Nucleotide
        :param nucl: Nucleotide object
        :return: tuple, in order of found nucl first
        """
        (n1_5, n1_3) = self.n1.to_tuple
        (n2_5, n2_3) = self.n2.to_tuple

        if nucl == n1_5:
            return n1_5, n2_3
        elif nucl == n2_3:
            return n2_3, n1_5
        elif nucl == n2_5:
            return n2_5, n1_3
        else:
            return n1_3, n2_5

    def onmodule(self, module):
        """
        The crossover does not have any built-in directionality, thus we do not know which pair of nucleotides
            is on which side of the crossover (and on which strand)
        Helps discern which NuclPair is on which module
        :param module: the module object
        :return: the NuclPair that is on module
        """
        m1 = self.n1.up()
        m2 = self.n2.up()
        if m1 == module:
            return self.n1
        elif m2 == module:
            return self.n2
        else:
            raise ValueError("This crossover does not exist on Ring({},{}).".format(module.bp, module.height))

    def __eq__(self, other):
        """
        Checks equivalence by matching all 4 nucleotide crossover points
        :param other: XoverPair
        :return: True/False
        """
        if type(other) != type(self):
            return False
        x1n1n5 = self.n1.n5
        x1n1n3 = self.n1.n3
        x1n2n5 = self.n2.n5
        x1n2n3 = self.n2.n3
        x1 = [x1n1n5, x1n1n3, x1n2n5, x1n2n3]
        x2n1n5 = other.n1.n5
        x2n1n3 = other.n1.n3
        x2n2n5 = other.n2.n5
        x2n2n3 = other.n2.n3
        x2 = [x2n1n5, x2n1n3, x2n2n5, x2n2n3]
        if all(nid in x2 for nid in x1) and all(nid in x1 for nid in x2):
            return True
        else:
            return False

    def get_nuclpair(self, nucl):
        """
        Gets the NuclPair nucleotides corresponding to the input Nucleotide
        :param nucl: Nucleotide object
        :return: tuple
        """
        (n1_5, n1_3) = self.n1.to_tuple
        (n2_5, n2_3) = self.n2.to_tuple

        if nucl == n1_5:
            return n1_5, n1_3
        elif nucl == n1_3:
            return n1_3, n1_5
        elif nucl == n2_3:
            return n2_3, n2_5
        else:
            return n2_5, n2_3

    def get_nucls(self):
        return self.n1.n5, self.n1.n3, self.n2.n5, self.n2.n3

    def __contains__(self, value):
        """
        Checks if Nucleotide value is present in the XoverPair
        :param value:
        :return:
        """
        if value in (self.n1.n5, self.n1.n3, self.n2.n5, self.n2.n3):
            return True
        else:
            return False

    def spans(self, m1, m2):
        """
        Returns True/False whether crossover spans the input modules
        :param m1: modules.Shape
        :param m2: modules.Shape
        :return: True/False
        """
        if issubclass(type(m1), modules.Shape) and issubclass(type(m2), modules.Shape):
            nt1 = self.n1.n5.get_top_module()
            nt2 = self.n2.n5.get_top_module()

        elif type(m1) == modules.JoinedStrand and type(m2) == modules.JoinedStrand:
            nt1 = self.n1.n5.get_top_strand()
            nt2 = self.n2.n5.get_top_strand()

        # Checks equality by their topology value (should be correctly set)
        if all(x in [m1, m2] for x in [nt1, nt2]):
            return True
        else:
            return False

    def __repr__(self):
        return "\n{:^15} {:^15} {}\n{:^15} {:^15}\n{:^15} {:^15} {} Bond Length: {}".format(
            self.n1.n5.numid, self.n1.n3.numid, self.n1.n3.get_top_module(), "|", "|", self.n2.n3.numid, self.n2.n5.numid, self.n2.n5.get_top_module(), round(self.bondlen(), 2))

    def bondlen(self):
        return np.linalg.norm(self.n2.coords - self.n1.coords)


    to_tuple = property(fget, fset)


class NuclPairPos:
    """
    Returns the angular position of a NuclPair
    ONLY FOR RING MODULES
    """
    def __get__(self, instance, owner):
        return mymath.getmdpt(mymath.nucl2angle(instance.n5), mymath.nucl2angle(instance.n3))
    
    def __set__(self, instance, value):
        raise AttributeError("Position cannot be set.")


class NuclPairTheta:
    """
    Returns the average theta angle value of two adjacent nucleotides of a NuclPair
    """
    def __get__(self, instance, owner):
        return mymath.getmdpt(instance.n5.theta, instance.n3.theta)
    
    def __set__(self, instance, value):
        raise AttributeError("Theta cannot be set.")


class NuclPairCoords:
    """
    Returns the midpoint of the centers of nucleotides in a NuclPair
    """
    def __get__(self, instance, owner):
        return (instance.n5.center+instance.n3.center)/2

    def __set__(self, instance, value):
        raise AttributeError("Coordinates cannot be set.")


class NuclPairBisectingVector:
    """
    Returns the bisecting vector of two nucleotides representing their average direction
    """
    def __get__(self, instance, owner):
        a = instance.n5.nucl_vec
        b = instance.n3.nucl_vec
        c = np.linalg.norm(a)*b + np.linalg.norm(b)*a
        return c

    def __set__(self, instance, value):
        raise AttributeError("BisectingVector cannot be set.")


class NuclPair:
    """
    Facilitate the format of a nucleotide pair as type tuple
    Define: a nucleotide pair is a pair of nucleotides that are adjacent
      on the same strand
    A nucleotide pair should be formatted:
    (<5' nucleotide>, <3' nucleotide>)

    2/17/2020: Removed condition that nucleotides have to be on same Module
                Modified to being on same Plane
    """
    pos = NuclPairPos()
    theta = NuclPairTheta()
    coords = NuclPairCoords()
    vector = NuclPairBisectingVector()

    def __init__(self, *nucl, nicked=False):
        self.nicked = nicked
        if self.chkformat(*nucl):
            (self.n5, self.n3) = nucl
            self.top = self.n5.up().up()  # to ringmodule
    
    def fget(self):
        """
        Gets NuclTuple
        :return: 5' and 3' Nucleotide objects
        """
        return self.n5, self.n3
    
    def fset(self, *value):
        """
        Sets NuclTuple
        :param value: 5' and 3' Nucleotide arguments
        :return: None
        """
        if self.chkformat(self, *value):
            (self.n5, self.n3) = value
            
    def nget(self):
        return self.n5.numid, self.n3.numid
    
    def nset(self):
        pass
    
    # returns ring
    def up(self):
        return self.top
        
    def chkformat(self, *value):
        """
        Checks if the input are Nucleotide objects
        :param value: Nucleotide objects
        :param nicked: True/False for adjacency/connectivity
        :return:
        """
        # Check that input is a 2-tuple
        try:
            (nucl5, nucl3) = value
        except:
            raise ValueError("{} does not match (5',3') nucleotide format.".format(value))

        # Check that both elements are Nucleotide objects
        if not (isinstance(nucl5, nucleotide.Nucleotide) and isinstance(nucl3, nucleotide.Nucleotide)):
            raise TypeError("Cannot set {} or {} to type 'Nucleotide'".format(type(nucl3), type(nucl5)))

        # Check for adjacency or connectivity
        if self.nicked:
            if nucl5.__strand3__ != nucl3 and nucl3.__strand5__ != nucl5:
                raise RuntimeError("Nucleotides must be adjacent to form a NuclPair.")
        else:
            if nucl5.toThree != nucl3 and nucl3.toFive != nucl5:
                raise RuntimeError("Nucleotides must be connected to form a NuclPair.")

        # Check that they are on the same strand
        if nucl5.get_top_strand() != nucl3.get_top_strand():
            raise AttributeError("Nucleotides do not belong to same strand. n5 [{}] n3 [{}]".format(
                nucl5.up().up(),
                nucl3.up().up()))

        return True
    
    to_tuple = property(fget, fset)
    numid = property(nget, nset)


class GapNuclPair(NuclPair):
    """
    Inherits from NuclPair for crossovers spanning gaps (half crossover placed on each side)
    """
    def __init__(self, *nucl):
        if self.chkformat(*nucl):
            (self.n5, self.n3) = nucl

        # # Break off the excess helix
        # # Scaffold
        # dnaconnector.nuclbreak(self.n5)
        # dnaconnector.nuclbreak(self.n3.__strand5__)
        # # Staple
        # dnaconnector.nuclbreak(self.n5.Comp.__strand5__)
        # dnaconnector.nuclbreak(self.n3.Comp)
        # # The strand's original absolute connections also need to be terminated
        # # Scaffold
        # dnaconnector.truebreak(self.n5)
        # dnaconnector.truebreak(self.n3.__strand5__)
        # # Staple
        # dnaconnector.truebreak(self.n5.Comp.__strand5__)
        # dnaconnector.truebreak(self.n3.Comp)

        # Emulate a cycle by spanning the gap for the staple and scaffold
        dnaconnector.trueconnect(self.n5, self.n3)
        dnaconnector.trueconnect(self.n3.Comp, self.n5.Comp)
        # Reconnect scaffolds (they will be disconnected later anyways to make the crossover)
        dnaconnector.nuclconnect(self.n5, self.n3)
        # Do not reconnect staples
        # TODO: Get rid of this, this is horrible
        dnaconnector.nuclconnect(self.n3.Comp, self.n5.Comp)

        # Removed constraint for gaps
        # if self.n5.toThree != self.n3 and self.n3.toFive != self.n5:
        #     raise RuntimeError("Nucleotides must be connected to form a NuclPair.")

        if self.n5.get_top_strand() != self.n3.get_top_strand():
            raise AttributeError("Nucleotides do not belong to same strand. n5 [{}] n3 [{}]".format(
                self.n5.up().up(),
                self.n3.up().up()))
        else:
            self.top = self.n5.up().up()  # to ringmodule

    def chkformat(self, *value):
        """
        Checks if the input are Nucleotide objects
        :param value:
        :return:
        """
        try:
            (nucl5, nucl3) = value
        except:
            raise ValueError("{} does not match (5',3') nucleotide format.".format(value))
        if not (isinstance(nucl5, nucleotide.Nucleotide) and isinstance(nucl3, nucleotide.Nucleotide)):
            raise TypeError("Cannot set {} or {} to type 'Nucleotide'".format(type(nucl3), type(nucl5)))
        # Removed constraint for gaps
        # if nucl5.toThree != nucl3 or nucl3.toFive != nucl5:
        #     raise ValueError("Nucleotides {} and {} are not connected or adjacent.".format(nucl5.numid, nucl3.numid))
        # else:
        #     return True
        return True


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
        :param xoverset: list of crossover.XoverPair
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


def filtervalid(xovers1, xovers2, **kwargs):
    """
    From a pair of lists of NuclPairs indicating valid crossover positions
    find any pair that has acceptable separation (< VALIDXOVERTHRESHBP)

    :param xovers1: list of NuclPair
    :param xovers2: list of NuclPair
    :return: List of XoverPair objects
    """
    def project_plane(x, angle):
        """
        :param x: ference crossover
        :param angle: angular difference between planes 1 and 2
        :return:
        """
        plane = x.up().up()
        # Project plane of crossover x onto flat plane
        xp_diff = x.coords - plane.origin
        if angle > 90:  # If the orientation is flipped, do not rotate
            xp_diff = np.concatenate((xp_diff[:1], rmatrix.rotate(xp_diff[1:], 0)), axis=0)
        else:  # Small rotations project to flat plane
            xp_diff = np.concatenate((xp_diff[:1], rmatrix.rotate(xp_diff[1:], -plane.yaw)), axis=0)
        return xp_diff

    def plane_spacing(x1, x2):
        """
        Determines the spacing applied between non-parallel modules
        :param x1: reference crossover 1
        :param x2: reference crossover 2
        :return: interhelical spacing if modules are projected to be parallel
        """
        if x1.up().bp != x2.up().bp:  # If different radius
            r1 = mymath.bp2radius(x1.up().bp)
            r2 = mymath.bp2radius(x2.up().bp)
            diff = abs(r1 - r2)
            return np.sqrt(np.power(config.INTERHELICAL, 2) - np.power(diff, 2))
        else:  # If same radius
            return config.INTERHELICAL

    threshold = kwargs.get("thresh", config.VALIDXOVERTHRESHBP)
    valid = []

    if not xovers1 or not xovers2:
        # raise NoValidXoverThreshold("Potential alignment positions are empty.")
        log.log_warning("Potential alignment positions are empty.")
        return valid
    for x1 in xovers1:  # match all pairs of potential crossover positions (NuclPairs)
        for x2 in xovers2:
            plane1 = x1.up().up()
            plane2 = x2.up().up()
            if plane1.yaw != plane2.yaw:
                '''Case: if shape axis curves
                NOT IMPLEMENTED FOR AXIALLY SYMMETRIC
                '''
                SPACING = plane_spacing(x1, x2)  # Set the spacing if planes are projected to be parallel
                # Project x1 and x2 modules each onto a flat plane
                yaw_angle = abs(mytrig.angle_diff(plane2.yaw, plane1.yaw))
                qx1c = project_plane(x1, yaw_angle)
                qx2c = project_plane(x2, yaw_angle)
                # After the projections, the modules will be overlapping so we need to manually space them one
                #   interhelical unit apart
                qx2c = qx2c + np.array([0, 0, SPACING])  # Space them one interhelical unit apart
                # Calculate the distance and filter
                dist = np.linalg.norm(qx2c - qx1c)
                if dist < threshold:
                    if all(n.strand_type() == 'stap' for n in [x1.n5, x1.n3, x2.n5, x2.n3]):
                        if all(not strandnav.neargap(n, config.GAP_SPACING_XOVER) for n in [x1.n5, x1.n3, x2.n5, x2.n3]):
                            valid.append(XoverPair(x1, x2))
                    elif all(n.strand_type() == 'scaf' for n in [x1.n5, x1.n3, x2.n5, x2.n3]):
                        valid.append(XoverPair(x1, x2))
                    else:
                        raise RuntimeError("Nucleotides are on an unknown strand.")
            else:
                '''Case: if shape axis is straight'''
                dist = np.linalg.norm(x2.coords - x1.coords)
                if dist < threshold:
                    if all(n.strand_type() == 'stap' for n in [x1.n5, x1.n3, x2.n5, x2.n3]):
                        if all(not strandnav.neargap(n, config.GAP_SPACING_XOVER) for n in [x1.n5, x1.n3, x2.n5, x2.n3]):
                            valid.append(XoverPair(x1, x2))
                    elif all(n.strand_type() == 'scaf' for n in [x1.n5, x1.n3, x2.n5, x2.n3]):
                        valid.append(XoverPair(x1, x2))
                    else:
                        raise RuntimeError("Nucleotides are on an unknown strand.")

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
    bondlen = kwargs.get("threshold", config.VALIDXOVERTHRESHBP) + config.THRESH_ADD
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
