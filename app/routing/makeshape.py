# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 19:18:40 2018

@author: Dan

4/2/2018:
    - makeshape.py create shapes as bidirectional graphs
    - nodes are the center positions
    - edges determine which nodes are connected
    - first finds coordinates, then places point objects at those coordinates
    
4/5/2018:
    - create freeform object
    - this generalizes the process of tracing a contour along points
    - encompasses line object
4/18/2018:
    - Remove bp2radius and univector math functions
6/21/2018:
    - Made theta private (removed self) from Circle3d
7/12/2018:
    - Added offtheta to initialization parameters of Circle3d for starting graph at specific angle
8/20/2018:
    - Minor change: change mymath label to mm
    - Update for consistency (Line), graph should be a list of coordinates
    - linetohelixbundle name changed to lineto6hb
    
"""
import math
import numpy as np
from .helper import (
    mymath,
    rmatrix as rm)
from app import config


def findclosestpoints(shape1, shape2):
    """
    Finds a pair of closest points on the line graphs of two shape objects to help determine the relative angle
        of the two shapes. Assumes that the two shapes are parallel.
        Since each shape has some granularity, search all pairs of points
    :param shape1: a Circle3d, Line, or Arc shape object
    :param shape2: a Circle3d, Line, or Arc shape object
    :return: tuple of two Point1d objects
    """
    # Initial value, pick a middle point to avoid edge cases
    point1 = shape1.graph[int(len(shape1.graph)/2)]
    point2 = shape2.graph[int(len(shape2.graph)/2)]
    mindist = np.linalg.norm(point2.coordinates - point1.coordinates)
    minpoints = (point1, point2)
    # Compare every pair of points but ignore edge cases within range 5
    for point1 in shape1.graph[5:-5]:
        for point2 in shape2.graph[5:-5]:
            distance = np.linalg.norm(point2.coordinates - point1.coordinates)
            if distance < mindist:
                mindist = distance
                minpoints = (point1, point2)
    return minpoints


class Point1d(object):
    """
    Describes a point in 3d
    Attributes define node location and connections
    """
    def __init__(self, next_point, prev_point, coordinates):
        """
        Initializes point with degree 2 connections and position
        :param next_point: point object
        :param prev_point: point object
        :param coordinates: 3d numpy array
        """
        self.next_point = next_point
        self.prev_point = prev_point
        self.coordinates = coordinates

        # error flags if node was not linked
        self.scaf = 'ERROR'
        self.stap = 'ERROR'

    def top(self, strand):
        """
        Returns the scaf or stap linked nucleotide
        """
        if strand == 'scaf':
            return self.scaf
        else:
            return self.stap

    def findnn(self, target_graph):
        """
        Finds the closest point to this point on the target module
        :param target_graph: Graph of taget module
        :return: Point1d on that module
        """
        point1 = self
        # Initial values
        point2 = target_graph[0]
        mindist = np.linalg.norm(point2.coordinates - point1.coordinates)
        minpoint = point2
        # Compare point1 to every point in the target graph
        for point2 in target_graph:
            distance = np.linalg.norm(point2.coordinates - point1.coordinates)
            if distance < mindist:
                mindist = distance
                minpoint = point2
        return minpoint


'''
symmetrical
'''


class Circle3d(object):
    """
    7/16/2018: add offtheta parameter to theta enumeration

    Creates a circle in 3D cartesian. Assumes it lays flat in xy-plane
    """
    def __init__(self, num_points, z, offtheta):
        self.num_points = num_points
        self.radius = mymath.bp2radius(num_points)
        self.normal = np.array([0, 0, 1])  # flat
        
        theta = []
        # populate all angles
        inc = 360.0/float(self.num_points)
        for n in range(self.num_points):
            theta.append((offtheta+n*inc) % 360.0)
        
        self.graph = []
        # populate all coordinate points
        for t in theta:
            x = self.radius*math.cos(math.radians(t))
            y = self.radius*math.sin(math.radians(t))
            # z already assigned
            # initialize unconnected nodes at coordinates into graph
            self.graph.append(Point1d(-1, -1, np.array([x, y, z])))
            
        # define explicit connection (instead of using indices)
        for i, node in enumerate(self.graph, start=0):
            node.next_point = self.graph[(i+1) % num_points]
            node.prev_point = self.graph[(i-1) % num_points]


class CircleOnPlane:
    """
    Independent function for now to create a circle on a defined plane rather than fixed Z
    TODO: Merge into Circle3d
    """
    def __init__(self, num_points, center, basis, offtheta=0):
        self.num_points = num_points
        self.radius = mymath.bp2radius(num_points)
        self.normal = np.array(basis[2])  # Equal to surface normal of plane

        theta = []
        # Populate all angle
        inc = 360.0 / float(self.num_points)
        for n in range(self.num_points):
            theta.append((offtheta + n * inc) % 360.0)

        self.graph = []
        # populate all coordinate points
        for t in theta:
            x = self.radius * math.cos(math.radians(t))
            y = self.radius * math.sin(math.radians(t))
            vx = x * basis[0]
            vy = y * basis[1]
            [x, y, z] = vx + vy + center
            # initialize unconnected nodes at coordinates into graph
            self.graph.append(Point1d(-1, -1, np.array([x, y, z])))

        # define explicit connection (instead of using indices)
        for i, node in enumerate(self.graph, start=0):
            node.next_point = self.graph[(i + 1) % num_points]
            node.prev_point = self.graph[(i - 1) % num_points]


'''
asymmetrical
'''


class Line(object):
    """
    Creates a line in 3D cartesian. Assumes it lays flat in xy-plane
    """
    # start, dir_vec must be vectors of type numpy.array
    def __init__(self, num_points, start_point, dir_vec, increment=config.AXIAL_RISE):

        # Error if not flat
        if dir_vec[2] != 0:
            raise ValueError("Only flat modules (z=0) are currently supported.")
        self.num_points = num_points
        self.normal = np.array([0, 0, 1])  # flat
        # convert to unit vector for good measure
        self.dir_vec = mymath.unitvector(dir_vec)

        self.start_point = start_point

        # inc = config.AXIAL_RISE  # axial rise
        inc = increment
        self.end_point = self.start_point + self.dir_vec * (num_points-1) * inc

        # make the points
        self.graph = []
        for n in range(num_points):
            self.graph.append(Point1d(-1, -1, start_point + n * config.AXIAL_RISE * self.dir_vec))

        # connect the points
        for n1, n2 in zip(self.graph[0:-1], self.graph[1:]):
            n1.next_point = n2
            n2.prev_point = n1


class LineTo6hb(object):
    # form new graphs in expected 6 helix bundle positions
    # put one line in, get one line out, otherwise routing is different
    # graph must stay list of points graph[points]
    # but output of 6 will be list of graphs with lists of points bundle[graph[points]]
    def __init__(self, vertices, strand_num):
        self.vertices = vertices
        self.strand_num = strand_num
        self.graph = []
        self.normal = np.array([0, 0, 1])  # flat
        strand_vec = [
            np.array([-1.3, 2.252]),
            np.array([-2.6, 0]),
            np.array([-1.3, -2.252]),
            np.array([1.3, -2.252]),
            np.array([2.6, 0]),
            np.array([1.3, 2.252])]
        for v in vertices:
            if v.next_point == -1:
                log.system("Warning: Hit a break!")
                break
            dir_vec = mymath.unitvector(v.next_point.coordinates - v.coordinates)
        
            bY = np.array([0, 0, 1])
            bX = np.cross(dir_vec, bY)
            
            self.graph.append(Point1d(-1, -1,
                                      strand_vec[self.strand_num][0] * bX +
                                      strand_vec[self.strand_num][1] * bY +
                                      v.coordinates))
        num_points = len(self.graph)
        # define explicit connection (instead of using indices)
        for i, node in enumerate(self.graph, start=0):
            node.next_point = self.graph[(i+1) % num_points]
            node.prev_point = self.graph[(i-1) % num_points]


class Arc:
    """Forms an arc shape"""
    def __init__(self, center, arcbp, init_angle, arc_angle, forced_circumference=None):
        """
        Initialize the arc
        :param center: x, y, z array
        :param arcbp: int length
        :param init_angle: float angle where the arc starts
        :param arc_angle: float angle range of the arc
        """
        self.init_angle = init_angle
        self.arc_angle = arc_angle
        self.num_points = arcbp
        section = arc_angle/360
        circumference = arcbp/section
        if forced_circumference:
            self.radius = mymath.bp2radius(forced_circumference)
        else:
            self.radius = mymath.bp2radius(circumference)
        self.normal = np.array([0, 0, 1])  # flat
        # arc length is 0.332
        # arc radius is self.radius
        # arc angle in radians is
        # aL = r*theta
        # theta = aL/r
        # populate all angles
        theta = np.linspace(init_angle, init_angle + arc_angle, self.num_points)

        self.graph = []
        # populate all coordinate points
        for t in theta:
            x = self.radius * math.cos(math.radians(t % 360))
            y = self.radius * math.sin(math.radians(t % 360))
            z = 0  # scale for center position later
            # initialize unconnected nodes at coordinates into graph
            self.graph.append(Point1d(-1, -1, np.array([x, y, z])))

        # shift to new center
        for node in self.graph:
            node.coordinates += center
        self.center = np.array(center)

        # connect the points
        for n1, n2 in zip(self.graph[0:-1], self.graph[1:]):
            n1.next_point = n2
            n2.prev_point = n1

    def get_endpoints(self):
        """
        Returns a point and direction vector for each endpoint of the arc
        Calculate the endpoint direction vector as an iteration of the last 2 vectors
            For points -3, -2, and -1,
                Calculate the angle difference between v(-3, -2) and v(-2, -1)
                Then rotate v(-2, -1) by that angle to get the endpoint direction vector
        :return: ((end1point, end1vector), (end2point, end2vector)), each is numpy.array
        """
        end1 = self.graph[0].coordinates
        end2 = self.graph[-1].coordinates
        # Find the last two adjacent vectors
        end1_v32 = self.graph[1].coordinates - self.graph[2].coordinates
        end1_v21 = self.graph[0].coordinates - self.graph[1].coordinates
        # Take only in-plane 2D x-y
        end1_v32 = end1_v32[:-1]
        end1_v21 = end1_v21[:-1]
        # Find the rotation angle
        a = np.arccos(np.dot(end1_v32, end1_v21)/(np.linalg.norm(end1_v32)*np.linalg.norm(end1_v21)))
        # Extrapolate the endpoint vector
        end1v = rm.rotate(end1_v21, a)
        end1v = np.array([end1v[0], end1v[1], 0])
        # Pad the end point
        end1 = end1 + end1v * config.AXIAL_RISE * 4

        # Find the last two adjacent vectors
        end2_v32 = self.graph[-2].coordinates - self.graph[-3].coordinates
        end2_v21 = self.graph[-1].coordinates - self.graph[-2].coordinates
        # Take only in-plane 2D x-y
        end2_v32 = end2_v32[:-1]
        end2_v21 = end2_v21[:-1]
        # Find the rotation angle
        a = np.arccos(np.dot(end2_v32, end2_v21) / (np.linalg.norm(end2_v32) * np.linalg.norm(end2_v21)))
        # Extrapolate the endpoint vector
        end2v = rm.rotate(end2_v21, a)
        end2v = np.array([end2v[0], end2v[1], 0])
        # Pad the end point
        end2 = end2 + end2v * config.AXIAL_RISE * 3

        return (end1, end1v), (end2, end2v)


class ArcOnPlane:
    """Forms an arc shape"""
    def __init__(self, center, arcbp, init_angle, arc_angle, origin, basis, forced_circumference=None):
        """
        Initialize the arc
        :param center: x, y, z array
        :param arcbp: int length
        :param init_angle: float angle where the arc starts
        :param arc_angle: float angle range of the arc
        """
        self.init_angle = init_angle
        self.arc_angle = arc_angle
        self.num_points = arcbp
        section = arc_angle/360
        if forced_circumference:
            forced_circumference = forced_circumference / section
            self.radius = mymath.bp2radius(forced_circumference)
        else:
            circumference = arcbp / section
            self.radius = mymath.bp2radius(circumference)
        self.normal = basis[2]
        # arc length is 0.332
        # arc radius is self.radius
        # arc angle in radians is
        # aL = r*theta
        # theta = aL/r
        # populate all angles
        theta = np.linspace(init_angle, init_angle + arc_angle, self.num_points)

        self.graph = []
        # populate all coordinate points
        for t in theta:
            x = self.radius * math.cos(math.radians(t % 360))
            y = self.radius * math.sin(math.radians(t % 360))
            vx = x * basis[0]
            vy = y * basis[1]
            [x, y, z] = vx + vy
            # initialize unconnected nodes at coordinates into graph
            self.graph.append(Point1d(-1, -1, np.array([x, y, z])))
        center = np.append(center, [0]) + origin
        # shift to new center
        for node in self.graph:
            node.coordinates += center
        self.center = np.array(center)

        # connect the points
        for n1, n2 in zip(self.graph[0:-1], self.graph[1:]):
            n1.next_point = n2
            n2.prev_point = n1

    def get_endpoints(self):
        """
        Returns a point and direction vector for each endpoint of the arc
        Calculate the endpoint direction vector as an iteration of the last 2 vectors
            For points -3, -2, and -1,
                Calculate the angle difference between v(-3, -2) and v(-2, -1)
                Then rotate v(-2, -1) by that angle to get the endpoint direction vector
        :return: ((end1point, end1vector), (end2point, end2vector)), each is numpy.array
        """
        end1 = self.graph[0].coordinates
        end2 = self.graph[-1].coordinates
        # Find the last two adjacent vectors
        end1_v32 = self.graph[1].coordinates - self.graph[2].coordinates
        end1_v21 = self.graph[0].coordinates - self.graph[1].coordinates
        # Take only in-plane 2D x-y
        end1_v32 = end1_v32[:-1]
        end1_v21 = end1_v21[:-1]
        # Find the rotation angle
        a = np.arccos(np.dot(end1_v32, end1_v21)/(np.linalg.norm(end1_v32)*np.linalg.norm(end1_v21)))
        # Extrapolate the endpoint vector
        end1v = rm.rotate(end1_v21, a)
        end1v = np.array([end1v[0], end1v[1], 0])  # WARNING: This is still fixed z
        # Pad the end point
        end1 = end1 + end1v * config.AXIAL_RISE * 4

        # Find the last two adjacent vectors
        end2_v32 = self.graph[-2].coordinates - self.graph[-3].coordinates
        end2_v21 = self.graph[-1].coordinates - self.graph[-2].coordinates
        # Take only in-plane 2D x-y
        end2_v32 = end2_v32[:-1]
        end2_v21 = end2_v21[:-1]
        # Find the rotation angle
        a = np.arccos(np.dot(end2_v32, end2_v21) / (np.linalg.norm(end2_v32) * np.linalg.norm(end2_v21)))
        # Extrapolate the endpoint vector
        end2v = rm.rotate(end2_v21, a)
        end2v = np.array([end2v[0], end2v[1], 0])  # WARNING: This is still fixed z
        # Pad the end point
        end2 = end2 + end2v * config.AXIAL_RISE * 3

        return (end1, end1v), (end2, end2v)



# freeform shape
class Freeform(object):
    def __init__(self, vertices):
        self.vertices = vertices
#        log.system(self.vertices)
        self.graph = []
        self.normal = np.array([0,0,1]) # flat
        inc = config.AXIAL_RISE # axial rise
        buffer_dist = 0
        # add the very first point
        self.graph.append(Point1d(-1, -1, vertices[0]['from']))
        for pair in vertices:
            dir_vec = pair['to']-pair['from']
            udir_vec = mymath.unitvector(dir_vec)
            max_dist = np.linalg.norm(dir_vec)
#            log.system("Going from",pair['from'],"to",pair['to'],"|Distance:",np.linalg.norm(max_dist))
            didx = 0 # distance index
            while True:
                if didx*inc+buffer_dist == 0: # catch error: haven't moved
                    didx+=1
                # the current point
                p = pair['from']+(didx*inc+buffer_dist)*udir_vec
                p_dist = np.linalg.norm(p-pair['from']) # distance traveled

                # make sure it is valid (not running off the endpoint edge)
                if p_dist>max_dist:
                    buffer_dist = -np.linalg.norm(max_dist-p_dist)
                    break
                
                self.graph.append(Point1d(-1, -1, p))
                didx+=1
        num_points = len(self.graph)
        # define explicit connection (instead of using indices)
        for i,node in enumerate(self.graph,start=0):
            node.next_point = self.graph[(i+1)%num_points]
            node.prev_point = self.graph[(i-1)%num_points]
        
        # check for cycle
        if not (self.vertices[0]['from'] == self.vertices[-1]['to']).all():
            self.graph[0].prev_point = -1
            self.graph[-1].next_point = -1

if __name__ == "__main__":
    bp = 128
    radius = mymath.bp2radius(bp)
    circ = Circle3d(bp, 0)
    
    bp = int(10.5*10)
    start_point = np.array([0, 0, 0])
    dir_vec = np.array([0, 1, 0])
    line = Line(bp, start_point, dir_vec)
    
    vertices = [{'from':np.array([0, 0, 0]),'to':np.array([10, 0, 0])}]
    ff = Freeform(vertices)
    newff = []
    for snum in range(6):
        newff.append(LineTo6hb(ff.graph, snum))