# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 18:33:00 2019

@author: Dan
"""

import numpy as np
import matplotlib

matplotlib.use('Agg')
from matplotlib import pyplot as plt
import shapely.geometry as geometry


# stl2ptcloud re-samples an STL mesh file
# input stl vertices (v)
# output ptsCloud which is 3D coordinate of each point
# Typically want around 100,000 points, but not much more than that, otherwise starts slowing down
def stl2ptcloud(v, MINPOINTS):
    import math

    # generate point cloud
    # print('Generating point cloud')
    ptsCloud = [];

    # current number of faces and points, each face has 3 points
    numFaces = len(v)
    numPoints = numFaces * 3

    # defaults
    sampledensity = 0
    stepsize = 1
    printstep = 100

    # downsample if too many points
    if numPoints > MINPOINTS:
        stepsize = math.floor(numPoints / MINPOINTS)
    # upsample if not enough points
    elif numPoints <= MINPOINTS:
        sampledensity = math.ceil(MINPOINTS / numPoints)

    for i in range(0, numFaces, stepsize):
        printIf = 0;
        if i % printstep == 0:
            # print('\tFace {} of {}:'.format(round(i / stepsize), round(len(v) / stepsize)))
            printIf = 1;
        [ptsCloud.append(_pts) for _pts in triPtsSample(v[i], sampledensity, printIf)]

    return ptsCloud


# TriPtsSample
# Converts triangulated face to sample points using Barycentric method
# http://blackpawn.com/texts/pointinpoly/
#   Input is ...
#       va,vb,vc 3 vertices of triangle
#       [old] density number of points to sample per area
#       number of additional random samples.
#   Output ...
#       each instance generates at least the 3 corner cases
#       generates at least one loop of 3 additional random points
def triPtsSample(vertices, samples, printIf):
    import random
    # save corner cases
    pts = [c for c in vertices]

    # define vectors
    (va, vb, vc) = vertices
    vecBA = vb - va
    vecCA = vc - va

    # number of sample points as function of area
    area = 0.5 * np.linalg.norm(np.cross(vecBA, vecCA))

    if printIf == 1:
        pass
        # print('\t\tSampling {} points in area {}'.format(samples, area))

    # iterate for desired number of sample points
    # each loop creates 3 sample points
    for i in range(samples):
        u = 1
        v = 1
        # 0<u,v or u+v<1 are required to not exceed triangle bounds
        while u + v > 1:
            u = random.random()
            v = random.random()
        # populate area
        pt = va + u * vecBA + v * vecCA
        pts.append(pt)

        # populate edge cases
        pt = va + u * vecBA
        pts.append(pt)
        pt = va + v * vecCA
        pts.append(pt)
    return pts


# For debugging the concave hull algorithm
# from descartes import PolygonPatch
# from matplotlib import pylab as pl
# def plot_polygon(polygon):
#     fig = pl.figure(figsize=(10,10))
#     ax = fig.add_subplot(111)
#     margin = .3
#
#     x_min, y_min, x_max, y_max = polygon.bounds
#
#     ax.set_xlim([x_min-margin, x_max+margin])
#     ax.set_ylim([y_min-margin, y_max+margin])
#     patch = PolygonPatch(polygon, fc='#999999', ec='#000000', fill=True, zorder=-1)
#     ax.add_patch(patch)
#     pl.savefig("test1")
#     return fig

# concave hull algorithm
# see: https://gist.github.com/dwyerk/10561690
def alpha_shape(points, alpha):
    from shapely.ops import cascaded_union, polygonize
    from scipy.spatial.qhull import Delaunay
    """
    Compute the alpha shape (concave hull) of a set
    of points.
    @param points: Iterable container of points.
    @param alpha: alpha value to influence the
        gooeyness of the border. Smaller numbers
        don't fall inward as much as larger numbers.
        Too large, and you lose everything!
    """
    if len(points) < 4:
        # When you have a triangle, there is no sense
        # in computing an alpha shape.
        return geometry.MultiPoint(list(points)).convex_hull

    coords = points
    tri = Delaunay(coords)
    triangles = coords[tri.vertices]
    a = ((triangles[:, 0, 0] - triangles[:, 1, 0]) ** 2 + (triangles[:, 0, 1] - triangles[:, 1, 1]) ** 2) ** 0.5
    b = ((triangles[:, 1, 0] - triangles[:, 2, 0]) ** 2 + (triangles[:, 1, 1] - triangles[:, 2, 1]) ** 2) ** 0.5
    c = ((triangles[:, 2, 0] - triangles[:, 0, 0]) ** 2 + (triangles[:, 2, 1] - triangles[:, 0, 1]) ** 2) ** 0.5
    s = (a + b + c) / 2.0
    areas = (s * (s - a) * (s - b) * (s - c)) ** 0.5
    circums = a * b * c / (4.0 * areas)
    filtered = triangles[circums < (1.0 / alpha)]
    edge1 = filtered[:, (0, 1)]
    edge2 = filtered[:, (1, 2)]
    edge3 = filtered[:, (2, 0)]
    edge_points = np.unique(np.concatenate((edge1, edge2, edge3)), axis=0).tolist()
    m = geometry.MultiLineString(edge_points)
    triangles = list(polygonize(m))
    # plot_polygon(cascaded_union(triangles))  # debugging, save png to local machine
    return cascaded_union(triangles), edge_points
