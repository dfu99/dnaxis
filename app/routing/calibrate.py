#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  8 20:23:38 2018

@author: dfu

calibrate.py
Methods for performing ring alignment

CHANGELOG:
    - 9/8/2018:
        - File created. Separated calibration functions from origami.py
    9/12/2018:
        Added top level functions from sym_origami.py
    9/13/2018:
        get_ring_alignment_difference
            did not account for wraparound in difference
    9/25/2018:
        make methods more explicitly specialized
        change ring alignment to calculate by xovers expected to be connected
    10/08/2018:
        deleted get_ring_alignment_difference
            obsolete: aligning by found crossover positions instead
        deleted calibrate_ring2
        removed top-level function.
            calibration now runs as call from route_staples
        
"""
import numpy as np
from . import (
        crossover)
from .helper import (
    log,
    mytrig,
    mymath)

from app import config


# xover_bp_pairs contains pairs of adjacent bases on the same strand
# [0] is the 5' base, [1] is the 3' base
# n - the reference angle position
# looks for the pair of bases that is closest in angle to n
# 10/08/2018: smallest difference between angles was not calculated correct
def findclosest(n, xover_bp_pairs):
    midpoint_positions = [mymath.get_midpoint(xp) for xp in xover_bp_pairs]
    positions = [abs(mytrig.angle_diff(mdpp, n)) for mdpp in midpoint_positions]
    index_min = np.argmin(positions)
    return xover_bp_pairs[index_min]


def matchtwist(nucl1, nucl2, **kwargs):
    """
    aligns the angle position of adjacent nucleotides at the endpoints of adjacent connected modules
    Uses nucl1 on strand 1 as reference to rotate nucl2 on strand 2
    :param nucl1: module 1 Nucleotide used as reference
    :param nucl2: module 2 Nucleotide on strand to calibrate
    """
    # Error handlers first
    # Check if basis vectors are opposite
    bX1 = nucl1.bX
    bX2 = nucl2.bX
    flagbasis = mymath.basis_dist(bX1, bX2, 1)
    if flagbasis:
        raise RuntimeError("The basis vectors {} and {} of adjacent modules at {} and {} is not consistent.".format(
            bX1, bX2, nucl1.numid, nucl2.numid))
    # Check if calibration is proceeding in 5'-3' direction
    # nucl1 of module1 should be the 5' base and nucl2 of module2 should be the 3' base
    if not (nucl1.toThree == -1 and nucl2.toFive == -1):
        log.system(nucl1.numid, nucl2.numid)
        log.system(nucl1.toThree, nucl2.toFive)
        try:
            log.system(nucl1.toThree.numid, nucl2.toFive.numid)
        except AttributeError:
            pass
        raise RuntimeError("Modules can only be connected in the 5'-3' direction (to avoid ambiguity), "
                           "or modules may already be connected.")
    # Otherwise proceed with the calibration by getting the angles of the adjacent nucleotides
    n1t = nucl1.theta
    n2t = nucl2.theta
    diff = mytrig.angle_diff(n1t, n2t)
    # If the transition is already within 30° (9 bp/t) to 40° (12 bp/t) angle difference, don't do anything
    # Direction is fixed, so the 3' angle should always lag the 5' angle
    #   and their difference should always be positive
    if 30 <= diff <= 40:
        m1apb = nucl1.up().angle_per_base
        m2apb = nucl2.up().angle_per_base
        avgrotation = (m1apb + m2apb) / 2
        log.debug(__name__, "Current: {} Optimal: {} but no twist calibration was performed.".format(diff, avgrotation))
        # print("Current: {} Optimal: {} but no twist calibration was performed.".format(diff, avgrotation))
        return
    # If the transition is not within an acceptable range, calibrate module2 using nucl2 as reference to the
    #   desired position
    else:
        # Set the new rotation as an average of the rotation per base pair of the two adjacent modules
        m1apb = nucl1.up().angle_per_base
        m2apb = nucl2.up().angle_per_base
        avgrotation = (m1apb+m2apb)/2
        # If a manual rotation has been set, override
        shift_angle = kwargs.get('angle', avgrotation)
        calibration_angle = (nucl1.theta - shift_angle) % 360
        # Logging
        log.debug(__name__, "Calibrating module to reference angle {}".format(calibration_angle))
        # print("Calibrating module to reference angle {}".format(calibration_angle))
        nucl2.up().twist(calibration_angle, nucl2)


def relax_twist(nucl1, nucl2, ref, span):
    """
    After modules are connected with matchtwist, there may be some remaining poor alignments.
    Use relax_twist to adjust the remaining misalignments to angle (ref) and distribute the excess to number
        of adjacent nucleotides (span) on each side of nucl1 and nucl2
    :param nucl1: Nucleotide of module 1
    :param nucl2: Nucleotide of module 2
    :return:
    """
    pass


# 20181229
# consolidated from calibrate_ring
def get_calibrated_ref(ring1, ring2):
    if config.CALIBRATIONMODE == "fixed":
        # fixed method
        # get the NuclPairs from the FullXover object that is on the correct ring
        _xovers1 = [x.onmodule(ring1) for x in ring1.objxovers]  # onmodule gets NuclPair
        # for a connection between ring n to n+1 in the path, use the
        # crossover from ring n-1 to n as a intermediary reference
        _ref1_xover = _xovers1[0]
        # find the position of the intermediary reference
        pos = mymath.get_midpoint(_ref1_xover)
        # interval of ring is how far an adjacent crossover should be
        # to maximize the spread of all crossovers
        ring2interval = 360/ring2.numxovers/2
        # move ring2interval away from pos
        refpos = (pos + ring2interval) % 360
        # get the valid crossover positions from ring1 to ring2
        xovers1 = ring1.getxoversto(ring2, 'stap')
        # and then find the closest valid crossover and use it as the reference
        ref1_xover = findclosest(refpos, xovers1)
    elif config.CALIBRATIONMODE == "dynamic":
        # dynamic method
        xovers1 = ring1.getxoversto(ring2, 'stap')
        ref1_xover = crossover.get_bestspace(xovers1)
    else:
        raise ValueError("CALIBRATIONMODE must be 'fixed' or 'dynamic'.")
    return ref1_xover


class Mixin:
    # =============================================================================
    # Calibration sub-functions
    # =============================================================================
    # calibrates ring1 and ring2 at a reference crossover, then returns that position
    def calibrate_ring(self, edge):
        ring1 = edge.directed['from']
        ring2 = edge.directed['to']
        (angle1to2, angle2to1) = mytrig.get_ringtoring_angle(ring1, ring2)
        # list is empty, no crossovers
        if not ring1.objxovers:
            # pick any valid crossover position from ring1 as the reference
            log.debug(__name__,
                      "No existing crossovers on Ring({},{}). Getting all possible crossovers for alignment.".format(
                          ring1.bp, ring1.height))
            xovers1 = ring1.getxoverstoring(ring2, 'stap')
            ref1_xover = xovers1[0]
        # list already has crossovers to the previous ring in the path
        else:
            log.debug(__name__,
                      "Found existing crossovers on Ring({},{}), using them as intermediary reference for alignment.".format(
                          ring1.bp, ring1.height))
            ref1_xover = get_calibrated_ref(ring1, ring2)

        # pick any valid crossover position from ring2 as the reference
        # this is ok because ring2 is assumed to be completely unconnected
        xovers2 = ring2.adjacent[ring1]['stxovers']
        ref2_xover = xovers2[0]
        # perform the alignment
        difference = self.get_xover_alignment_difference(ref1_xover, ref2_xover)
        ring2.rotate(difference)
        # return the aligned crossover position as reference
        xover_pair = crossover.XoverPair(ref1_xover, ref2_xover)
        return xover_pair

    def get_xover_alignment_difference(self, xover1, xover2):
        """
        Checks what the angular difference is between xover pairs
        xover1 and xover2 are NuclPair objects
        :param xover1: XoverPair(object) from ring1 is stationary
        :param xover2: XoverPair(object) from ring2 being rotated
        :return: rotation value in degrees (float)
        """
        x1a = mymath.get_midpoint(xover1)
        x2a = mymath.get_midpoint(xover2)
        diff = mytrig.angle_diff(x1a, x2a)
        return diff

    def clean_twist(self):
        """
        Looks through all modules for adjacent nucleotides that still have some misalignment
        :return:
        """
        for strand in self.get_all_strands():
            for nucl5 in strand:
                try:
                    nucl3 = nucl5.__strand3__
                    if not nucl5.up().up() is nucl3.up().up():  # not on same module, so connected
                        b5 = nucl5.bX
                        b3 = nucl3.bX
                        if mymath.basis_dist(b5, b3, 1):  # if different basis, normalize
                            nt2 = nucl5.theta
                            nt1 = mymath.lh2rh2(nucl3.theta)
                            diffnt = abs(mytrig.angle_diff(nt2, nt1))  # check the difference
                        else:
                            nt2 = nucl5.theta
                            nt1 = nucl3.theta
                            diffnt = abs(mytrig.angle_diff(nt2, nt1))
                        if not 30 <= diffnt <= 40:  # outside of acceptable angle_per_base
                            if diffnt > 40:
                                ref = 40
                            else:
                                ref = 30
                            relax_twist(nucl5, nucl3, ref, 10)
                        else:
                            raise RuntimeError("Something else went wrong.")
                except AttributeError:
                    # hit a nick, which there shouldn't be, since strand should still be cyclic
                    raise AttributeError("Twist relaxation hit a nick."
                                         "Should only be used while strands are still cyclic.")

