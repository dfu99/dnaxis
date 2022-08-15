#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 19:03:37 2018

@author: dfu

collide.py

description:
    Functions to be called by scaffold crossover routing functions (scxover.py).
    After a scaffold crossover is placed, it may be too close to a staple crossover.
    If it is within the same ring, the staple crossover should simply be removed.
    If it collides with a staple crossover in a different ring.
"""

from .helper import mymath, mytrig, strandnav

from operator import itemgetter


class Mixin:
    # top level
    # cleans up all staple xovers that collide with nearby scaffold xovers
    # perhaps unnecessary to consolidate
    def clean(self):
        pass # TODO
    
    # gets all scaffold crossovers
    def get_allscafxovers(self):
        pass # TODO
    
    # gets all staple crossovers
    def get_allstapxovers(self):
        pass # TODO
        
    # current version that checks only the edge between two adjacent rings
    # check if xoverpair is colliding with an existing xover
    def check_colliding(self, x):
        # get which ring each NuclPair of the Xover is on
        r1 = x.n1.up()
        r2 = x.n2.up()
        
        # get the already present crossovers on those rings
        r1x = r1.objxovers
        r2x = r2.objxovers
        r1r2x = r1x+r2x
        
        # filter for crossovers only on the same rings as Xover x
        rx = [_x for _x in r1r2x if xoveronrings(_x, r1, r2)]
        
        # find closest crossover on same rings
        cx = self.closest_fx(x, rx)
        
        # only undo it if it is too close - NOT ACTUALLY IMPLEMENTED
        # define too close as within one turn
        cxobj = cx[1]
        cxobj.undo()

    # find the closest existing xover to a potential xover
    def closest_fx(self, xp, fxlist):
        distxp2fx = [min((strandnav.xp_distfromxover(xp)), fx) for fx in fxlist]
        sortdist = sorted(distxp2fx, key=itemgetter(0))
        # log.out(__name__, sortdist)
        return sortdist[0]

    
# get average position of FullXover or XoverPair
def posxo(fx):
    return mymath.getmdpt(fx.n1.pos, fx.n2.pos)


# helper function to determine if an xover lies between two rings
def xoveronrings(xover, r1, r2):
    n1 = xover.n1
    n2 = xover.n2
    if nponrings(n1, r1, r2) and nponrings(n2, r1, r2):
        return True
    return False


# helper function to determine if a nuclpair is on either of two adjacent rings
def nponrings(np, r1, r2):
    if np.up() == r1 or np.up() == r2:
        return True
    return False