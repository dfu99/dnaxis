# -*- coding: utf-8 -*-
"""
Created on Sun Oct 14 21:02:22 2018

@author: Dan

symmesh.py translates MATLAB code from fSymMesh.m to Python
to remove dependency on MATLAB engine import which only seems to run on Python 3.5

CHANGELOG:
    10/14/2018
        File created
    8/30/2019
        Rewrite MATLAB code into Python
"""
import numpy as np
from stl import mesh
import os
import math
from shapely.geometry import LineString, MultiPoint
from app.routing.helper import mymath
from . import helixgeom as hxg
from . import mypoly as myp
from more_itertools import pairwise


def readstl(filepath):
    # Read the CAD file
    #    print("Opening {}".format(filename))
    v = mesh.Mesh.from_file(os.path.join(filepath))
    # print(v.points)
    # f is an array. each element is a 3-element vector designating the
    # vertices of a face
    return v


def stl2ptcloud(v):
    # Resample to Point cloud
    # use vertices of STL to populate more points
    pts = myp.stl2ptcloud(v.vectors, 10000)
    # remove duplicate rows
    # not sure why there are duplicates in the first place
    pts = np.unique(pts, axis=0)
    # Save a 3D figure here
    return np.array(pts)


def flatten(pts, center='y'):
    # Flattened 2D point cloud
    if center == 'y':
        x = pts[:, 0]
        y = pts[:, 1]
    elif center == 'z':
        x = pts[:, 1]
        y = pts[:, 2]
    else:
        x = pts[:, 2]
        y = pts[:, 0]
    # concatenate into n x 2 vector array
    return np.c_[x, y]


def fitboundary(pts, alpha):
    x = pts[:, 0]
    y = pts[:, 1]
    # ADJUST THE ALPHA UNTIL GETTING A GOOD FIT
    concave_hull, edge_points = myp.alpha_shape(pts, alpha=alpha)

    fullpts = [np.array([p[0], p[1]]) for p in list(concave_hull.exterior.coords)]

    for p in pairwise(fullpts):
        p0 = p[0]
        p1 = p[1]
        if p0[0] < 0 <= p1[0]:
            break
        else:
            fullpts = fullpts[1:] + [fullpts[0]]

    # cut in half
    halfpts = np.array([np.array([p[0], p[1]]) for p in fullpts if p[0] >= 0])
    fullpts = np.array(fullpts)

    (hx, hy) = (halfpts[:, 0], halfpts[:, 1])
    (fx, fy) = (fullpts[:, 0], fullpts[:, 1])
    halfsec = LineString(np.c_[hx, hy])
    fullsec = LineString(np.c_[fx, fy])

    return halfsec, fullsec


def openshape(xsec):
    # openshape removes flat, horizontal edges
    # Samples 'flatres' points, then takes them two at a time
    # If their vector difference has an angle that is at most flatdeg away from 0 or 180 degrees
    # They are considered horizontal and skipped.
    # @parameter, input, xsec
    #   LineString, traces half the shape
    # @output
    #   LineString with flat edges removed

    # Number of points we'll sample for corners
    flatres = 360
    # Threshold below which we consider a horizontal line
    flatdeg = 1
    # Split the line into flatres even segments
    delflat = MultiPoint([xsec.interpolate((i / (flatres - 1)), normalized=True) for i in range(0, flatres)])
    dnodes = np.array([np.array(pt.coords[0]) for pt in list(delflat)])

    # Ignore each flat edge
    dpts = []
    for d in pairwise(dnodes):
        v = d[1] - d[0]
        pol = mymath.cart2pol(v[0], v[1])
        sinmag = abs(math.sin(pol[1]))
        if sinmag < math.sin(flatdeg):
            pass
        else:
            dpts.append(d[0])
            dpts.append(d[1])

    dpts = np.array(dpts)

    return LineString(dpts)


def fithelices(user_numrings, user_xoverCount, user_maxnt, xsec, use_tpx=True, MINCIRCUMFERENCE=72, interhelical=2.6,
               interhelicalmaxdist=3.0, MINTPX=2, MAXTPX=5, MINBPT=9, MAXBPT=12):
    while True:
        # Divide into ring helix positions
        # lineSplit

        splitter = MultiPoint(
            [xsec.interpolate((i / (user_numrings - 1)), normalized=True) for i in range(0, user_numrings)])
        nodes = np.array([np.array(pt.coords[0]) for pt in list(splitter)])

        avgdist = np.mean(mymath.distanceBtwnNodes(nodes))
        # Scale to interhelical distance
        scalingFactor = interhelical / avgdist
        nodes = hxg.scaleshape(nodes, scalingFactor)

        # print("Distance Check:\n")
        for i in range(1, len(nodes)):
            currPos = nodes[i]
            lastPos = nodes[i - 1]
            distance = np.linalg.norm(currPos - lastPos)
            # print("{}".format(distance))

        # Convert to ring objects
        # initialize
        rings = []
        tpx = MINTPX
        # print("Node creation:\n")
        for i in range(len(nodes)):
            x = nodes[i, 0]
            h = nodes[i, 1]
            # print(x, h)

        for i in range(len(nodes)):
            x = nodes[i, 0]
            h = nodes[i, 1]

            rings.append(hxg.Ring(mymath.radius2bp(x) + MINCIRCUMFERENCE, tpx, user_xoverCount, h))

        # Scale smallest ring to minimum accepted circumference
        # Calculate minimum accepted circumference
        if use_tpx:
            minimumRingCircumference = MINTPX * MINBPT * user_xoverCount
        else:
            minimumRingCircumference = MINCIRCUMFERENCE

        # Find the smallest ring
        smallestCircumference = rings[0].bp
        for i in range(user_numrings):
            if rings[i].bp < smallestCircumference:
                smallestCircumference = rings[i].bp

        # scale rings if necessary
        if smallestCircumference < minimumRingCircumference:
            rings = hxg.scalerings(rings, minimumRingCircumference - smallestCircumference)

        # Iterate until stable
        tpxrange = np.arange(MINTPX, MAXTPX + 1, 1)
        for r in rings:
            # print("Stabilizing ring {}".format(str(r)))
            # ring-specific
            stable = False
            while (stable == False):
                # if valid tpx can be applied, ring is stable
                # otherwise, unstable
                for t in tpxrange:
                    bpt = r.bbx / t
                    if bpt <= MAXBPT and bpt >= MINBPT:
                        r.resetRing(r.bp, t, r.numxovers, r.height)
                        stable = True
                        break
                    else:
                        stable = False
                # if unstable and bpt>12, double the crossovers
                if stable == False and r.bbx / MAXTPX > 12:
                    r.resetRing(r.bp, r.tbx, r.numxovers * 2, r.height)
                # if unstable, scale the ring by one unit of numx
                elif stable == False:
                    # print("I'm scaling!")
                    r.resetRing(r.bp + r.numxovers, r.tbx, r.numxovers, r.height)

        # Proof for interhelical distance violation
        for j in range(len(rings) - 1):
            R0 = rings[j].radius
            h0 = rings[j].height
            R = rings[j + 1].radius
            h = rings[j + 1].height

            dist = np.linalg.norm(np.array([R - R0, h - h0]))
            if dist > interhelicalmaxdist:
                pass
                # print("Warning! Helices distanced too far!")

        # Proof for scaffold length restrictions
        totalntlen = 0
        for r in rings:
            totalntlen += r.bp
        if totalntlen <= user_maxnt:
            # print("Current length {} is within acceptable {}. Continuing to output.".format(totalntlen, user_maxnt))
            break
        elif totalntlen > user_maxnt:
            # print("Current length {} is not within {}. Reducing number of rings and trying again.".format(totalntlen,
            #                                                                                               user_maxnt))
            user_numrings -= 1
    return rings
