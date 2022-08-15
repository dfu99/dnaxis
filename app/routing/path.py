#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  8 20:15:27 2018

@author: dfu
path.py

CHANGELOG:
    - 9/8/2018:
        - File created. Separated path functions from origami.py
    9/12/2018:
        - Added top-level functions from sym_origami.py
    9/27/2018:
        Added Edge class
            Added __init__
            Added method to_set
            Added method equals
    9/28/2018:
        Added method path_to_edges
"""
from numpy import linalg
import numpy as np
import itertools
from .helper import mytrig, mymath, log
from . import makeshape, modules


class Edge:
    """# 20180927
    Edge class is the graph representation of adjacent modules
    """

    def __init__(self, module1, module2):
        self.vertices = (module1, module2)
        self.directed = {'from': module1, 'to': module2}
        self.thresh = self.get_thresh()

    def __eq__(self, other):
        """
        Checks if another edge's vertices are the same as its own
        :param other: Edge(object)
        :return: TRUE if same, FALSE if not
        """
        if type(other) != type(self):
            return False
        elif all(v in self.vertices for v in other.vertices):
            return True
        else:
            return False

    def __hash__(self):
        return hash(repr(self))

    def __get__(self, instance, owner):
        return list(self.vertices)

    def toString(self):
        module1 = self.directed['from']
        module2 = self.directed['to']
        s = "from module {} to module {}".format(module1, module2)
        return s

    def __contains__(self, item):
        if item in self.vertices:
            return True
        else:
            return False

    def __repr__(self):
        m1 = self.vertices[0]
        m2 = self.vertices[1]
        return m1.note+"<->"+m2.note

    def get_thresh(self):
        """
        Returns either shortest or average distance to adjacent module depending on whether they are parallel or not
        :return: int value
        """
        module1 = self.vertices[0]
        module2 = self.vertices[1]
        plane1 = module1.get_top_plane()
        plane2 = module2.get_top_plane()
        if plane1.yaw != plane2.yaw:
            return module1.shortest_dist_to(module2)
        else:
            return module1.avg_dist_to(module2)


class EdgeSet:
    """
    Custom methods for sets of Edge objects
    Uses undirected edges
    """

    def __init__(self, li):
        """
        :param li: Input should be a list of Edge objects
        """
        self.value = mymath.to_set(li)

    def __iter__(self):
        return iter(list(self.value))

    # adds an Edge to the EdgeSet if it is not already there
    def add(self, edge):
        if self.value:  # If existing elements, avoid adding a duplicate
            for existing_edge in self.value:
                if existing_edge == edge:
                    return
        self.value.append(edge)

    # checks for membership
    def contains(self, edge):
        for v in list(self.value):
            if v == edge:
                return True
        return False

    def getsorted(self):
        """
        sorts the edges such that we go from inside - outside, bottom - top
        :return: list of sorted Edges
        """
        srtlst = []
        tmp = list(self.value)
        while tmp:
            currmin = tmp[0]
            for t in tmp:
                # if it's lower, takes precedence
                con1 = t.directed['from'].height < currmin.directed['from'].height
                # if it's on the same plane, but inside
                con2 = t.directed['from'].height == currmin.directed['from'].height
                con3 = t.directed['from'].bp < currmin.directed['from'].bp
                if con1:  # if lower
                    currmin = t
                elif con2 and con3:  # if on the same plane, inside
                    currmin = t
            srtlst.append(currmin)
            tmp.remove(currmin)
        return srtlst

    # remove looks for edge in set
    # looks for complement if KeyError
    def remove(self, edge):
        try:
            self.value.remove(edge)
        except KeyError:
            for comp_edge in self.value:
                if comp_edge.equals(edge):
                    self.value.remove(comp_edge)
                    return
            raise KeyError("Path {} does not exist in {}.".format(edge.toString(), self.__class__))
        finally:
            log.debug(__name__, "Removed path {} from the edge set.".format(edge.toString()))


class Mixin:
    # =============================================================================
    # Top-level structural functions
    # =============================================================================
    # assumes planes are already sorted in grid
    # then connections should only be by adjacent nodes
    # each node has degree 4 edges
    def eval_connections(self):
        log.out(__name__, "Parsing planes for connected modules.")
        grid = self.grid
        gridpairs = mymath.gridpairs(grid)
        for pair in gridpairs:
            ((r1i, r1j), (r2i, r2j)) = pair
            ring1 = self.grid[r1i][r1j]
            ring2 = self.grid[r2i][r2j]
            (angle1to2, angle2to1) = mytrig.get_ringtoring_angle(ring1, ring2)
            ring1.setadjacent(ring2, angle1to2)
            ring2.setadjacent(ring1, angle2to1)

    def auto_connections(self, threshold):
        """
        Automatically generates the dict of connected modules
        :return: connections dict
        """
        connections = {}
        all_modules = self.get_modules()
        for module in all_modules:  # For all E=(v1, v2), populate all v1 into dict
            connections[self.get_module_index(module)] = []
        for m1, m2 in itertools.combinations(all_modules, 2):  # Pair up all modules
            (p1, p2) = makeshape.findclosestpoints(m1.shape, m2.shape)
            if linalg.norm(p1.coordinates - p2.coordinates) < threshold:
                mid1 = self.get_module_index(m1)
                mid2 = self.get_module_index(m2)
                connections[mid1].append(mid2)
                connections[mid2].append(mid1)
        return connections

    def prune_connections(self, connections, pct_overlap_thresh, dist_thresh=3.0):
        """
        For every pair of connected modules, check that is there is enough overlap to make crossovers
        :param connections: connections dict
        :param pct_overlap_thresh: Threshold percentage of points sufficiently overlapping
        :param dist_thresh: Threshold distance between nodes
        :return: connections dict with insufficiently overlapping pairs removed
        """
        if pct_overlap_thresh <= 0:
            return connections
        elif pct_overlap_thresh > 1:
            log.log_warning("Bad overlap percentage submitted to __name__.")
            return connections
        # Iterate through every pair of connections between module1 and module2
        for key in connections:
            module1 = self.get_module_by_index(key)
            for idx in connections[key]:
                module2 = self.get_module_by_index(idx)
                # Initial value, pick a middle point to avoid edge cases
                point1 = module1.shape.graph[int(len(module1.shape.graph) / 2)]
                point2 = module2.shape.graph[int(len(module2.shape.graph) / 2)]
                mindist = np.linalg.norm(point2.coordinates - point1.coordinates)
                # Track distance to closest point on module2 for every point on module1
                arr_mindist = []
                for point1 in module1.shape.graph:
                    for point2 in module2.shape.graph:
                        distance = np.linalg.norm(point2.coordinates - point1.coordinates)
                        if distance < mindist:
                            mindist = distance
                    arr_mindist.append(mindist)
                # Threshold and check how many points are within range of each other
                pct_overlap = [1 if x <= dist_thresh else 0 for x in arr_mindist]
                pct_overlap = sum(pct_overlap)/len(pct_overlap)
                if pct_overlap < pct_overlap_thresh:
                    log.system("Overlap {}% between modules {}{} and {}{}. Should be removed.".format(pct_overlap * 100, str(module1), key, str(module2), idx))
                    connections[key].remove(idx)
                    connections[idx].remove(key)
                else:
                    log.system("Overlap {}% between modules {}{} and {}{}. Acceptable.".format(pct_overlap * 100, str(module1), key, str(module2), idx))
        return connections

    def connect_path_multi(self):
        log.out(__name__, "Evaluating connection pathway.")
        sequence = []  # Initialize output variable
        this_ring = self.planes[0].modules[0]  # pick the first location
        self._weaveout = True  # initialize for beginning on the inside
        while True:
            sequence.append(this_ring)
            adjacency_list = this_ring.adjacent
            # remove visited indicies
            next_list = [x for x in adjacency_list if x not in sequence]

            # if all indices were removed, we're at the terminus
            # else, figure out what to do next
            if not next_list:  # next_list is empty == False
                return sequence
            else:
                this_ring = self.route_next_ring(this_ring, next_list)

    def apply_connections(self, connections):
        for c in connections:
            ring1 = self.planes[c - 1].modules[0]
            for r in connections[c]:
                ring2 = self.planes[r - 1].modules[0]
                (angle1to2, angle2to1) = mytrig.get_ringtoring_angle(ring1, ring2)
                ring1.setadjacent(ring2, angle1to2)
                ring2.setadjacent(ring1, angle2to1)
        return

    def apply_asym_connections(self, connections):
        """
        Temporary function. To be merged into apply_connections once all formats are normalized again.
        Biggest difference is that apply_connections assumes there is only one module per plane.
        Thus connections input is only one value denoting the plane index.
        The new connections format should be a Dict, with Keys as the (plane, module) index tuple and the
            the Values is a list of connecting modules also in that format
        :param connections: Dict
        :return: None
        """
        for key in connections:
            try:
                module1 = self.planes[key[0]].modules[key[1]]
            except IndexError:
                raise IndexError("Could not find {}".format(key))
            for val in connections[key]:
                try:
                    module2 = self.planes[val[0]].modules[val[1]]
                except IndexError:
                    raise IndexError("Could not find {}".format(val))
                # This method of calculating angles may CHANGE along the path of the shape
                # (p1, p2) = makeshape.findclosestpoints(module1.shape, module2.shape)
                # (x, y, z) = p2.coordinates - p1.coordinates
                # (rho, theta, phi) = mymath.cart2sph(x, y, z)
                # phi = math.degrees(phi)

                plane1 = module1.up()
                plane2 = module2.up()
                # Only for same diameter RingModules, pitch = 0, yaws different
                if plane1.yaw != plane2.yaw and \
                        plane1.pitch == plane2.pitch == 0 and \
                        type(module1) == type(module2) == modules.RingModule and \
                        module1.bp == module2.bp:
                    # print("DEBUG.apply_asym_connections: Confirm yaw being run. Code 0xE32")
                    # Convert yaw to a unit vector
                    yaw1 = plane1.yaw
                    (v1x, v1y) = (np.cos(yaw1), np.sin(yaw1))
                    v1 = mymath.unitvector(np.array([v1x, v1y]))
                    yaw2 = plane2.yaw
                    (v2x, v2y) = (np.cos(yaw2), np.sin(yaw2))
                    v2 = mymath.unitvector(np.array([v2x, v2y]))
                    # Find the midpoint
                    yaw0 = (yaw1 + yaw2) / 2
                    (v0x, v0y) = (np.cos(yaw0), np.sin(yaw0))
                    v0 = mymath.unitvector(np.array([v0x, v0y]))

                    # Calculate an offset
                    offset = np.degrees(np.arccos(np.dot(v2, v0)))
                    if yaw2 > yaw1:
                        o1 = offset
                        o2 = -offset
                    else:
                        o1 = -offset
                        o2 = offset

                    r1 = module1.radius
                    r2 = module2.radius
                    z1 = 0
                    z2 = 2.6
                    d = linalg.norm(np.array([r1, z1]) - np.array([r2, z2]))
                    phi = mytrig.get_angle(r1, r2, z1, z2, d)
                    module1.setadjacent(module2, phi+o1)

                else:  # For parallel modules, no yaw
                    # print("DEBUG.apply_asym_connections: This should not have run. Code 0xE31")
                    try:  # For ArcModule and RingModule
                        r1 = module1.radius
                        r2 = module2.radius
                        z1 = module1.height
                        z2 = module2.height
                        d = linalg.norm(np.array([r1, z1]) - np.array([r2, z2]))
                        phi = mytrig.get_angle(r1, r2, z1, z2, d)
                    except AttributeError:  # For LineModule
                        (p1, p2) = makeshape.findclosestpoints(module1.shape, module2.shape)
                        (_, _, phi) = mymath.cart2sph(*(p2.coordinates - p1.coordinates))
                        phi = np.degrees(phi)

                    # Goes through angles of nucleotide pair of the helix and checks if phi intersects it
                    module1.setadjacent(module2, phi)

    def apply_pathway(self, pathway):
        """
        Implementation currently for symmetric structures only
        :param pathway: List of integers (1-indexed) corresponding to order of planes (0-indexed)
        :return:
        """
        sequence = [self.planes[r - 1].modules[0] for r in pathway]
        # # re-normalize all the directions
        # # this doesn't do anything because the directions are already set when the rings are made
        # # Obsolete
        # _dirBit = True
        # for ring in sequence:
        #     ring.dirBit = _dirBit
        #     _dirBit = not _dirBit
        return sequence

    def apply_asym_pathway_edges(self, pathway):
        """
        Temporary method for testing methods specific to the asymmetric example 'Clover'.
            To be consolidated into apply_pathway later
        Directly initiates the Edge graph pathway since the nodes are not consecutive
        :param pathway: List of 2-Tuple denoting an Edge (module 1, module 2)
                Each module is denoted with another 2-Tuple (plane index, module index)
        :return: List of Edge objects
        """
        edges_list = []
        for edge in pathway:
            mdx1 = edge[0]
            module1 = self.planes[mdx1[0]].modules[mdx1[1]]
            if edge[1]:
                mdx2 = edge[1]
                module2 = self.planes[mdx2[0]].modules[mdx2[1]]
                edges_list.append(Edge(module1, module2))
        return edges_list

    # converts a pathway formatted as nodes
    # to edge format
    def path_to_edges(self):
        edges_list = []
        for module1, module2 in zip(self.pathway[0:-1], self.pathway[1:]):
            edges_list.append(Edge(module1, module2))
        return edges_list

    # =============================================================================
    # Path routing sub-functions
    # =============================================================================
    # go through several cases of determining which adjacent module to route to next
    # looking only at visited indices, prefer in plane
    # then check the nearest out of plane        
    def route_next_ring(self, this_ring, next_list):
        if any(x.height == this_ring.height for x in next_list):
            # find the closest ring in plane
            in_plane = [x for x in next_list if x.height == this_ring.height]
            next_ring = get_closest_ring(this_ring, in_plane)
            if next_ring.bp > this_ring.bp:  # determine weaving direction
                self._weaveout = True
            else:
                self._weaveout = False
            this_loc = [this_ring.bp, this_ring.height]
            log.debug(__name__, f"(Route in) @ %s. Potential connections: %s" %
                      (str(this_loc), str([[x.bp, x.height] for x in in_plane])))
        else:
            next_ring = get_border_ring(next_list, self._weaveout)
            this_loc = [this_ring.bp, this_ring.height]
            log.debug(__name__, f"(Route out) @ %s. Potential connections: %s" %
                      (str(this_loc), str([[x.bp, x.height] for x in next_list])))
        return next_ring

    def partition_from_removed(self, removed):
        """
        Separate the modules in the origami based on where scaffold crossover were removed for separated scaffolds
        :param removed: List of removed edges
        :return: List of lists of modules
        """
        partitions = []  # List of lists of modules
        this_part = []  # List of modules
        for edge in self.pathway:
            mdx1 = edge[0]
            module1 = self.planes[mdx1[0]].modules[mdx1[1]]
            mdx2 = edge[1]
            module2 = self.planes[mdx2[0]].modules[mdx2[1]]
            this_edge = Edge(module1, module2)
            for module in module1.up().modules:
                if module not in this_part:
                    this_part.append(module)
            if removed.contains(this_edge):
                partitions.append(this_part)
                this_part = []
            else:
                for module in module2.up().modules:
                    this_part.append(module)
        partitions.append(this_part)
        return partitions


def get_closest_ring(ring1, ringlist):
    """
    Returns the ring from list based on closest radius
    """
    sorted_ringlist = sorted(ringlist, key=lambda x: abs(x.bp - ring1.bp), reverse=False)
    return sorted_ringlist[0]


# returns the largest or smallest ring from list depending on side arg
def get_border_ring(self, plane, side):
    try:
        sorted_plane = sorted(plane, key=lambda x: x.bp, reverse=side)
        return sorted_plane[0]
    except:
        log.out(__name__, "get_border_ring argument side expects binary.")


def edge_list_to_path_list(li):
    """
    Converts a list of source and target vertices to a list of ordered vertices
    :param li: List of 2-tuple (source, target vertex). Each vertex is a 2-tuple (plane id, module id)
    :return: Collapsed ordered list
    """
    new_list = []
    # Add each source vertex
    for tup in li[:-1]:
        source = tup[0]
        new_list.append(source)
    # Add the source and target of the last element
    source = li[-1][0]
    target = li[-1][1]
    new_list.append(source)
    new_list.append(target)

    return new_list
