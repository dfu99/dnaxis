#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 16:39:58 2018

@author: dfu

Last modified by Daniel Fu (Reif Lab, Duke University) on 4/5/2020

Name        : mymath
Project     : cadaxisdna
Description : Some basic calculations and actions
Interpreter : Python 3.7.4

Changelog:
    4/19/2018:
        File created
        Some basic functions added
    9/10/2018:
        File organization, moved mymath functions from sym_origami.py to here
    9/12/2018:
        Added nucldist
    9/13/2018:
        Adjusted get_midpoint to account for wraparound
    10/08/2018
        Consolidated get_basis_vec from DNAHelix
        Fixed math in nucl2angle
    10/10/2018
        Switch get_midpoint input to expect NuclPair
"""

import numpy as np
import math
import itertools
from . import mytrig
from app import config

"""
SPECIAL TYPES
"""


class Tape:
    """
    Queue where all elements are kept in the collection
    """

    def __init__(self, li):
        self.value = li
        self.reset()

    def reset(self):
        self._ind = -1

    def valueAt(self, v):
        return self.value[v]

    def forward(self):
        self._ind += 1
        try:
            return self.valueAt(self._ind)
        except IndexError:
            raise StopIteration

    def backward(self):
        self._ind -= 1
        # still indexerror if we've done something to go too far the other way
        return self.valueAt(self._ind)

    def current(self):
        return self.valueAt(self._ind)


class FakeSet:
    """
    Ordered Set. Removes the randomness of implementing collections as Set, but still removes duplicate elements.
    """

    def __init__(self, l):
        if type(l) == list:
            self.value = l
        else:
            raise ValueError("{} is not a list".format(l))

    # only add if element is not already in list                
    def add(self, v):
        if v not in self.value:
            self.value.append(v)
        else:
            pass

    def __contains__(self, v):
        if v in self.value:
            return True
        else:
            return False

    def __len__(self):
        return len(self.value)

    def __sub__(self, other):
        return [a for a in self.value if a not in other]

    def __iter__(self):
        self.i = 0
        return self

    def __next__(self):
        if self.i < len(self.value):
            result = self.value[self.i]
            self.i += 1
            return result
        else:
            raise StopIteration

    def __getitem__(self, item):
        return self.value[item]

    def __setitem__(self, key, val):
        self.value[key] = val


"""
DNA ORIGAMI SPECIFIC MATH FUNCTIONS
"""


def bp2radius(c):
    """
    Converts circumference in bps to radius in nm
    :param c: Circumference in terms of number of base pairs
    :return: Returns the radius of the ring
    """
    return c * config.AXIAL_RISE / (2 * math.pi)


def radius2bp(r):
    """
    Converts radius in nm to circumference in bps
    :param r: Radius of the arc in nm
    :return: Circumference of the full circle
    """
    return 2 * np.pi * r / config.AXIAL_RISE


def nucl2angle(nucl):
    (x, y, z) = nucl.center
    if x == 0:
        if y < 0:
            angle = 270
        else:
            angle = 90
    elif y == 0:
        if x < 0:
            angle = 180
        if x > 0:
            angle = 0
    else:
        angle = math.degrees(math.atan(y / x))
        if x > 0:
            pass
        else:
            angle += 180
    return angle % 360


def nucldist(n1, n2):
    """
    returns Euclidean space distance between backbone molecules of two nucleotides
    :param n1: type Nucleotide
    :param n2: type Nucleotide
    :return: float distance
    """
    n1pos = np.array([n1.x, n1.y, n1.z])
    n2pos = np.array([n2.x, n2.y, n2.z])
    dist = np.linalg.norm(n1pos - n2pos)
    return dist


"""
GENERAL MATH FUNCTIONS
"""


def distanceBtwnNodes(nodes):
    numNodes = len(nodes)
    distances = np.zeros(numNodes-1)

    for i in range(1, numNodes):
        lastNode = nodes[i-1]
        currPos = nodes[i]
        distances[i-1] = np.linalg.norm(currPos - lastNode)
    return distances


def roundnearest(x, base=5):
    return base * round(x/base)


def get_precision(flt):
    """
    Gets the precision in decimal points of a number
    :param flt: float value
    :return: int decimal points
    """
    val = flt
    precision = 0
    if val % 1 == 0:  # Checks if integer
        return 0
    while val % 1 != 0:
        val = val * 10
        precision += 1
    return precision


def to_set(edges_list):
    """
    to_set takes in a list of values
    Checks for duplicates
    Returns as a set but uses list typing
    """
    edges_set = edges_list
    for e1, e2 in itertools.combinations(edges_list, 2):
        if e1 == e2:
            edges_set.remove(e2)
    return edges_set


def mean(li):
    """
    # calculate the average of a 1-D list
    :param li:
    :return:
    """
    return sum(li) / len(li)


def unitvector(v):
    """
    helper function that converts vector to its unit vector
    :param v: type numpy.array; any dimension vector
    :return: type numpy.array unit vector
    """
    return v / np.linalg.norm(v)


# expects NuclPair object
def get_midpoint(nuclpair):
    (bp1, bp2) = nuclpair.to_tuple
    p1 = nucl2angle(bp1)
    p2 = nucl2angle(bp2)
    return getmdpt(p1, p2)  # 20181026


def flipangle(angle):
    return (360 - angle) % 360


def basis_dist(b1, b2, dist):
    """
    Checks if two basis vectors are farther than dist
    :param b1:
    :param b2:
    :param dist:
    :return:
    """
    if np.linalg.norm(b2 - b1) > dist:
        return True
    return False


def getmdpt(a, b):
    """
    expects angle in degrees a, b. order does not matter.
    missing citation for where this algorithm was taken from
    :param a: angle in degrees
    :param b: angle in degrees
    :return: bisecting angle in degrees
    """
    p = max(a, b) % 360
    q = min(a, b) % 360

    if p - q > 180:
        d = p - 360
        a = (d + q) / 2
        if a < 0:
            a = a + 360
    else:
        a = (p + q) / 2
    return a


def unpack(np):
    """
    20181024 method created
    unpack a NuclPair object into angles
    :param np: type NuclPair
    :return: p1 is 5' nucleotide rotational position,
     p2 is 3' nucleotide rotational position"""
    (bp1, bp2) = np.to_tuple
    p1 = nucl2angle(bp1)
    p2 = nucl2angle(bp2)
    return p1, p2


def cycle_list(_list, first_index):
    """
    cycles the list such that the specified index is now first
    :param _list: list to cycle
    :param first_index: index to cycle to
    :return: cycled list
    """
    newlist = _list[first_index:] + _list[:first_index]
    return newlist


def get_basis_vec(shapenorm, dirnorm):
    """
    Find basis vectors for nucleotide plane normal to direction dir_vec
    :param shapenorm: Normal 3-d vector
    :param dirnorm: Direction 3-d vector
    :return:
    """

    bY = shapenorm  # Y in nucleotide plane
    bX = np.cross(dirnorm, bY)  # X in nucleotide plane
    bX = bX / np.linalg.norm(bX)
    return bX, bY


def lstlenmin(lst):
    """
    returns the shortest length element of a list
    :param lst: List of Lists
    :return: Shortest List
    """
    try:
        lmin = sorted(lst, key=lambda x: len(x))[0]
    except IndexError:
        raise IndexError("Input was empty.")
    return lmin


# assumes all items are unique
def index2d(lst, val):
    """
    Returns the row, column index looking for a value in a 2D list
    :param lst:
    :param val:
    :return:
    """
    for i, l in enumerate(lst):
        for j, v in enumerate(l):
            if v == val:
                return i, j
    raise ValueError("Value '{}' does not exist in 2d list.".format(val))


# returns the adjacent nodes in a grid array    
def gridpairs(lst):
    pairs = []
    for i in range(len(lst)):
        for j in range(len(lst[i])):
            try:
                lst[i][j + 1]
                pairs.append(((i, j), (i, j + 1)))
            except IndexError:
                pass
            try:
                lst[i + 1][j]
                pairs.append(((i, j), (i + 1, j)))
            except IndexError:
                pass
    return pairs


def distrolist(lst, num):
    """
    Provide num evenly spaced positions to insert into a list
    :param lst: List with at least num elements
    :return: List of evenly distributed indices
    """
    if len(lst) < num:
        raise ValueError("List is not long enough to be divided into {} positions.".format(num))
    m = len(lst)
    n = num
    return [i * m // n + m // (2 * n) for i in range(n)]


def roundeach(lst, dec):
    """
    Rounds each element in a list to the input number of decimals
    :param lst: list or tuple iterable of flats
    :param dec: int
    :return: list
    """
    out = []
    for each in lst:
        out.append(round(each, dec))
    return out


def min_to_bpt(r):
    """
    Calculates the turns and error to a given ring circumference to most closely adhere to nominal twist
    :param r: ring circumference
    :return: (turns, error)
    """
    best = {'dist': None, 'turns': None, 'i': None}
    for turns in range(1, round(r/10.6) * 2):
        for i in range(-5, 6, 1):
            dist = abs((r + i)/turns - 10.6)
            try:
                if dist < best['dist']:
                    best['dist'] = dist
                    best['turns'] = turns
                    best['i'] = i
            except TypeError:
                best['dist'] = dist
                best['turns'] = turns
                best['i'] = i
    return best


def smaller_module(module1, module2):
    """
    Returns smaller module
    :param module1: modules.Shape
    :param module2: modules.Shape
    :return: modules.Shape
    """
    if module1.bp <= module2.bp:
        return module1
    else:
        return module2


def base_to_module_angle(base, module):
    """
    Calculates the angle from base to target module
    :param base: Nucleotide
    :param module: modules.RingModule
        - Currently not supported for other modules.Shape child classes
    :return:
    """
    # Use the 5' node
    # But it shouldn't matter because all points on the graph should be more or less equidistant to another
    # point on the target graph.
    node = base.node
    # Find the closest point to node on target module
    target_node = node.findnn(module.shape.graph)
    # Calculate the center to center vector
    target_vector = unitvector(target_node.coordinates - node.coordinates)
    # Build the basis transition vector
    bx = base.bX
    by = base.bY
    bz = base.duvec
    b = np.c_[bx, by, bz]
    # Project the center-center vector onto the nucleotide's basis to the XY plane
    proj_target_vector = unitvector(project_3dvector_to_2dbasis(target_vector, bx, by))
    # Find the XY plane angle
    target_angle = custombasisangle(b, proj_target_vector)
    return target_angle


def base_to_module_dist(base, module):
    """
    Calculates the distance from base to target module
    :param base: nucleotide.Nucleotide
    :param module: modules.RingModule
        - Currently not supported for other modules.Shape child classes
    :return:
    """
    # Use the 5' node
    # But it shouldn't matter because all points on the graph should be more or less equidistant to another
    # point on the target graph.
    node = base.node
    # Find the closest point to node on target module
    target_node = node.findnn(module.shape.graph)
    # Calculate the center to center vector
    target_vector = target_node.coordinates - node.coordinates
    unit_target_vector = unitvector(target_vector)
    # Build the basis transition vector
    bx = base.bX
    by = base.bY
    bz = base.duvec
    b = np.c_[bx, by, bz]
    # Project the center-center vector onto the nucleotide's basis to the XY plane
    basis_target_vector = unitvector(project_3dvector_to_2dbasis(unit_target_vector, bx, by))
    # Then project the original target_vector onto vector in the basis to get its length
    # proj_a (b), b = target_vector; a = basis_target_vector
    proj_target_vector = projection(target_vector, basis_target_vector)
    # Find the magnitude of the vector
    target_magnitude = np.linalg.norm(proj_target_vector)
    return target_magnitude


def delete_chars(s, c_list):
    """
    Replaces multiple characters from a string with an empty string
    :param s: string
    :param c_list: characters
    :return: new string
    """
    for c in c_list:
        s = s.replace(c, "")
    return s


"""
Linear algebra
"""


def project_3dvector_to_2dbasis(v, bx, by):
    """
    Projects a 3d vector onto the plane defined by orthonormal basis vectors bx and by
    :param v: 3d vector
    :param bx: 3d vector
    :param by: 3d vector
    :return: (v1', v2'), the coefficients of the vector in the basis (bx, by)
    """
    # Project the vector on the basis vectors, thus into the plane
    vx = projection(v, bx)
    vy = projection(v, by)

    return vx + vy


def projection(b, a):
    """
    Project vector b onto a
    :param b: 3d vector
    :param a: 3d vector
    :return: 3d vector
    """
    return np.dot(a, b) / np.linalg.norm(a) ** 2 * a


def custombasisangle(bmatrix, vec):
    """
    Calculates the angle of vec if projected into a 2D plane defined by basis vectors in bmatrix
    :param bmatrix: 3x3 matrix
    :param vec: 3d Vector projected onto the plane
    :return: Returns an angle, in degrees
    """
    invbmatrix = np.linalg.inv(bmatrix)
    t = np.matmul(invbmatrix, vec.transpose())
    if round(t[2], 6) != 0:  # non-zero k coefficient means vec was not projected into the 2D plane
        raise ValueError("The input vector is not in the same 2D plane as defined by the basis. z={}.".format(t[2]))
    (x0, y0, x1, y1) = (0, 0, t[0], t[1])
    # print("DEBUG math domain error", x0, x1, y0, y1)
    theta = mytrig.get_angle(x0, x1, y0, y1, 1)
    return theta


def edgelist2nodelist(edgelist):
    """
    Converts a list of edges with initial and final node to a sequence of ordered nodes
    Must satisfy line path property, each node only has one outbound and inbound path
    :param edgelist:
    :return:
    """
    nodelist = [edgelist[0][0], edgelist[0][1]]
    for edge in edgelist[1:]:
        node1 = edge[0]
        node2 = edge[1]
        if node1 not in nodelist:
            nodelist.append(node1)
        if node2 not in nodelist:
            nodelist.append(node2)
    return nodelist


def linkedlisttoadjmtx(linkedlist, options='zigzag'):
    """
    Converts a linked list to an adjacency matrix

    :param linkedlist: defined as dict with key and list of adjacent objects
    :param options: 'zigzag' or 'layer' determines how the adjacency matrix is indexed, which will affected
        what direction the later pathway algorithm prioritizes
    :return: adjacency matrix
    """
    if options == 'layer':
        indices = tuple([k for k in linkedlist.keys()])
        print("Indices order:", indices)

    elif options == 'zigzag':
        indices = []
        indices.append(list(linkedlist.keys())[0])
        while True:
            print("Per iteration:", indices)
            planeidx = indices[-1][0]
            for j in linkedlist.keys():
                if j not in indices and j[0] == planeidx:
                    indices.append(j)
                    continue
                for k in linkedlist.keys():
                    if k not in indices:
                        indices.append(k)
                        continue
            break

    else:
        raise ProcessLookupError("Linked list to adjacency matrix conversion has no", options, "option.")
    dmxn = len(indices)
    am = [[0]*dmxn for i in range(dmxn)]
    for idx1, tup in enumerate(indices):
        for adj in linkedlist[tup]:
            idx2 = indices.index(adj)
            am[idx1][idx2] = 1
    return am


"""
COORDINATE SYSTEM CONVERSION
"""


def cart2cyl(x, y, z):
    rho = np.sqrt(x ** 2 + y ** 2)
    phi = np.arctan2(y, x)
    z = z
    return rho, phi, z


def cyl2cart(rho, phi, z):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    z = z
    return x, y, z


def cart2sph(x, y, z):
    """
    Converts 3d cartesian coordinates to spherical
    :param x: float
    :param y: float
    :param z: float
    :return: (rho, theta, phi) 3-tuple float
    """
    rho = math.sqrt(x ** 2 + y ** 2 + z ** 2)
    theta = np.arctan2(y, x)
    # if y >= 0:
    #     phi = np.arccos(z/rho)
    # else:
    #     phi = 2*np.pi-np.arccos(z/rho)
    phi = np.arctan2(np.sqrt(x ** 2 + y ** 2), z)
    if theta <= 0:
        phi = 2 * np.pi - phi

    # Conform to CADAxiSDNA directions
    # phi = lh2rh(phi)

    return rho, theta, phi


def lh2rh(phi):
    """
    Specific to CADAxiSDNA to normalize right-handed spherical coordinates to left-handed cartesian
    :param phi: float angle, radians
    :return: adjusted float angle, radians
    """
    # Custom transformtion for CADAxiSDNA to conform to left-handed angles
    phi = (abs(2 * math.pi - phi) + math.pi / 2) % (2 * math.pi)
    return phi


def lh2rh2(phi):
    """
    Specific to CADAxiSDNA to normalize right-handed spherical coordinates to left-handed cartesian
    :param phi: float angle, degrees
    :return: adjusted float angle, degrees
    """
    # Custom transformtion for CADAxiSDNA to conform to left-handed angles
    phi = (abs(360 - phi) + 180) % 360
    return phi


def array_cart2cyl(x, y, z):
    arr_rho = []
    arr_phi = []
    arr_h = []
    for i in range(len(x)):
        _x = x[i]
        _y = y[i]
        _z = z[i]
        (_r, _p, _h) = cart2cyl(_x, _y, _z)
        arr_rho.append(_r)
        arr_phi.append(_p)
        arr_h.append(_h)
    return arr_rho, arr_phi, arr_h


def array_cyl2cart(rho, phi, h):
    arr_x = []
    arr_y = []
    arr_z = []

    for i in range(len(rho)):
        _rho = rho[i]
        _phi = phi[i]
        _h = h[i]
        (_x, _y, _z) = cyl2cart(_rho, _phi, _h)
        arr_x.append(_x)
        arr_y.append(_y)
        arr_z.append(_z)
    return arr_x, arr_y, arr_z


def cart2pol(x, y):
    """Converts 2D cartesian coordinates into 2D polar
    :param x, y: float coordinates
    :return: tuple(rho, phi)
    """
    rho = np.sqrt(x ** 2 + y ** 2)
    phi = np.arctan2(y, x)
    return rho, phi


def pol2cart(rho, phi):
    """
    Converts 2D polar coordinates to 2D cartesian
    :param rho: magnitude
    :param phi: angle (in radians)
    :return: x, y
    """
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return x, y


def array_cart2pol(x, y):
    arr_rho = []
    arr_phi = []
    for i in range(len(x)):
        _x = x[i]
        _y = y[i]
        (_r, _p) = cart2pol(_x, _y)
        arr_rho.append(_r)
        arr_phi.append(_p)
    return np.array(arr_rho), np.array(arr_phi)


def array_pol2cart(rho, phi):
    arr_x = []
    arr_y = []

    for i in range(len(rho)):
        _rho = rho[i]
        _phi = phi[i]
        (_x, _y) = pol2cart(_rho, _phi)
        arr_x.append(_x)
        arr_y.append(_y)
    return np.array(arr_x), np.array(arr_y)


if __name__ == "__main__":
    # sample_list = [i for i in range(20)]
    # distribution = distrolist(sample_list, 5)
    # print(distribution)

    pass
