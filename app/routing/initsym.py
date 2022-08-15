# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 19:14:55 2018

@author: Dan

initsym.py
Split from initorigami.py for functions for only symmetrical structures

CHANGELOG:
    7/12/2018:
        - File split
        - draw_sym reduced to importing data, editing for layers, converting to ring object, the adding to origami
    8/21/2018:
        - Change routing to follow plane-by-plane separation
"""

from . import makeshape
from . import modules

from .filehandler import symfile

'''
### STRUCTURE INITIALIZATION
Calls shape creation methods
Then calls nucleotide populating methods
Creates origami object in designed shape
'''


def draw_sym(filename):
    raw_data = symfile.getfile(filename) # import as floats
    planes = init_planes(raw_data)
    link_planes(planes)
    return planes

# =============================================================================
# PLANE CREATION
# =============================================================================


def init_planes(data):
    """
    function initPlanes
    creates all planes where rings exist
    :param data: data is a 2D list of raw input data
    :return: a plane object with object modules added to it
    """
    # dirBit = True # track 5' or 3' direction
    nucl_idx_cnt = 1  # begin a count for nucleotides
    planes = []
    for line in data:
        # format parameters
        bp = int(line[0])
        tbx = int(line[1])
        xovers = int(line[2])
        height = float(line[3])
        dirBit = line[6]
        
        # make planes
        thisplane = modules.Plane(float(height))
        # make ring
        ring = modules.RingModule(bp, tbx, xovers, height, dirBit)
        # make circle shape
        shape_obj = makeshape.Circle3d(bp, height, 0)
        # apply shape
        ring.applyshape(shape_obj, nucl_idx_cnt)
        # set hierarchy
        ring.settop(thisplane)
        # increment nucleotide counter
        nucl_idx_cnt += ring.bp*2
        # add ring to plane
        thisplane.add_module(ring)
        # add plane to list
        planes.append(thisplane)
        # flip direction
        # dirBit = not dirBit
    return planes


# set connection data for planes
def link_planes(data):
    for i, plane in enumerate(data):
        if not i == 0:
            plane.adjacent.append(data[i-1])
        if not i == len(data)-1:
            plane.adjacent.append(data[i+1])


