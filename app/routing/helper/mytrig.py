# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 15:57:22 2018

@author: Dan
Changelog:
    4/1/2018:
        - angle wrapping function angle_between
    4/2/2018:
        - consolidate angle-related functions here in angle.py
        - write connection comprehension functions connection_comp
    5/7/2018:
        - Bug in angle_between:
        a1<=b and a2<=b should be
        a1<=b and b<=a2
        - Change the compared value from b -> n for less ambiguity
    5/10/2018:
        - write new function eval_xover_pos that checks nucleotides in the ring for valid crossover positions
        - rename connection_comp (connection comprehension) to eval_connections
            - also change layerdata to origami for consistency
    9/8/2018:
        - moved some functions from origami.py
    9/23/2018:
        - improvements to get_angle
            - added distance argument d
            - fixed hypotenuse is wrong
    10/08/2018:
        added angle_diff
"""

"""  
### ANGLE FUNCTIONS
"""
from numpy import linalg
import numpy as np
import math
from . import log


def get_angle(r0, r1, z0, z1, d):
    """
    Determine the angle between two adjacent rings.
    Only requires the circumference so it is independent of creating the routing.
    :param r0: Radius of initial ring (float)
    :param r1: Radius of target ring (float)
    :param z0: Height of initial ring (float)
    :param z1: Height of target ring (float)
    :param d: Distance between the rings, for the hypotenuse (float)
    :return: Angle from the initial ring to the target ring (float, unit degrees)
    """
    # find angle of r0 TO r1
    # r0 is FROM ring, r1 is TO ring
    # TODO: Add cases for z-axis position, 8/3/2017: assumes top to bottom
    h = d  # fixed distance to next ring (the hypotenuse)
    dR = abs(r1-r0)

    z0 = round(z0, 12)
    z1 = round(z1, 12)
    r0 = round(r0, 12)
    r1 = round(r1, 12)
    if z0 == z1:  # planar case
        if r0 < r1:
            return 0
        else:
            return 180
    elif z0 > z1:
        if r0 > r1:
            a = 90 + math.degrees(math.asin(dR/h))
        elif r0 < r1:
            a = 90 - math.degrees(math.asin(dR/h))
        else:
            a = 90  # for syntax purposes only
        return 360-a
    elif z0 < z1:
        if r0 < r1:
            a = 90 + math.degrees(math.asin(dR/h))
        elif r0 > r1:
            a = 90 - math.degrees(math.asin(dR/h))
        else:
            a = 90  # for syntax purposes only
        return 180-a
    else:
        return -1


def angle_between(a1, n, a2):
    """
    # determines if angle b lies between angles a1, a2
    # July 16 2018: Only works for acute angles
    # ref https://www.xarg.org/2010/06/is-an-angle-between-two-other-angles/
    # 20181024 added ValueError if a1, a2 is are colinear
    # 20181127 changed from <= and >= to <= and >
    :param a1:
    :param n:
    :param a2:
    :return:
    """
    if abs(angle_diff(a1, a2)) == 180:
        raise ValueError("Angles are colinear, whether {} is between {} and {} is ambiguous.".format(n, a1, a2))
    (a1, a2) = sort_acute(a1, a2)
    n = (360 + (n % 360)) % 360
    a1 = (360 + a1) % 360
    a2 = (360 + a2) % 360
    if a1 < a2:
        # print(n,"between",a1,"and",a2,"?",a1<=n and n<a2)
        return a1 <= n and n < a2
    # print(n,"between",a1,"and",a2,"?",a1<=n or n<a2)
    return a1 <= n or n < a2


def vector_between(a, b, c):
    """
    Checks if vector B lies between vectors A and C
    https://stackoverflow.com/questions/13640931/how-to-determine-if-a-vector-is-between-two-other-vectors
    Assume vectors lie in same plane
    :param a:
    :param b:
    :param c:
    :return:
    """
    return np.cross(a, b) * np.cross(a, c) >= 0 and np.cross(c, b) * np.cross(c, a) >= 0


def sort_acute(a1, a2):
    """
    Returns the order by which the positive angle from a1 to a2 is acute
    :param a1: angle, degrees
    :param a2: angle, degrees
    :return: tuple
    """
    t1 = (a1-a2) % 360
    t2 = (a2-a1) % 360
    if t1 > t2:
        return a1, a2
    else:
        return a2, a1


def acute_diff(a, b):
    """
    Finds the acute angle difference between angles a and b
    :param a: type float; unit degrees
    :param b: type float; unit degrees
    :return: d1 = a - b if d1<d2 or d2 = b - a if d2<d1
    """
    d1 = angle_diff(a, b)
    d2 = angle_diff(b, a)
    if d1 < d2:
        return d1
    return d2


def angle_diff(a, b):
    """
    Finds the acute angle difference between two angles, accounting for wraparound
    d = a (param1) - b (param2)
    :param a: type float; unit degrees
    :param b: type float; unit degrees
    :return: d = a - b; type float; unit degrees
    """
    d = a - b
    angle = (d+180) % 360-180
    return angle


def get_ringtoring_distance(ring1, ring2):
    """
    This is a simplification of a closest points algorithm, since shapes are assumed to be concentric circles
    Specific to CADAxisSDNA symmetric design
    :param ring1: Initial module.RingModule
    :param ring2: Target module.RingModule
    :return: float, Euclidean distance between these two ring modules
    """
    (r1, z1) = ring1.getposition_exact()
    (r2, z2) = ring2.getposition_exact()
    distance = linalg.norm(np.array([r1, z1])-np.array([r2, z2]))
    return distance


def get_ringtoring_angle(ring1, ring2):
    """
    This is a simplification of finding the angle between closest points on concentric circles
    Specific to CADAxiSDNA symmetric design
    :param ring1: Initial module.RingModule
    :param ring2: Target module.RingModule
    :return: float, degrees angle between these two ring modules
    """
    (r1, z1) = ring1.getposition_exact()
    (r2, z2) = ring2.getposition_exact()
    d = linalg.norm(np.array([r1, z1])-np.array([r2, z2]))
    angle1to2 = get_angle(r1, r2, z1, z2, d)
    angle2to1 = get_angle(r2, r1, z2, z1, d)
    return angle1to2, angle2to1


if __name__ == "__main__":
    import mymath
    # test cases
    r0 = mymath.bp2radius(96)
    r1 = mymath.bp2radius(144)
    z0 = 1.0
    z1 = 2.15
    # print(get_angle(r0, r1, z0, z1, 2.81))
