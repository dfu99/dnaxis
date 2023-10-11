"""
Functions for handling crossover discovery and placement that are specific to ASYMMETRIC mode designs
"""

from ..helper import mymath, rmatrix
from app import config
from app.routing.motifs import XoverPair
import numpy as np


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


"""
Case: if shape axis curves
                NOT IMPLEMENTED FOR AXIALLY SYMMETRIC
"""
def filtervalid(xovers1, xovers2, threshold):
    valid = []
    for x1 in xovers1:  # match all pairs of potential crossover positions (NuclPairs)
        for x2 in xovers2:
            plane1 = x1.up().up()
            plane2 = x2.up().up()
            if plane1.yaw != plane2.yaw:
                plane1 = x1.up().up()
                plane2 = x2.up().up()
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

    return valid
