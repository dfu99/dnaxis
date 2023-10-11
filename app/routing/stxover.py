#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stxover.py
Methods acting upon the staple crossovers
Created on Thu Oct 11 19:48:31 2018
As fork from crossover.py

@author: dfu

changelog:
    11/21/2018
        max(ring1.numxovers, ring2.numxovers) changed to min
"""

from . import crossover, heuxover
from app import config
from .helper import log
import math


class Mixin:
    """
    =============================================================================
    Staple stuff
    =============================================================================
    top level staple crossover function
    20180929: Changed number of xovers to be determined by
      the most xovers that can be added instead of ring1.xovers
    """
    def route_staples(self):
        log.out(__name__, "Adding staple crossovers.")
        # track edges for non-pathway edges
        unvisited_edges = self.get_edges()
        # iterate through pathway
        for edge in self.pathway_edges:
            ring1 = edge.directed['from']
            ring2 = edge.directed['to']
            log.out(__name__, "Connecting staple for Edge( Ring({},{}), Ring({},{}) ).".format(
                ring1.bp, ring1.height,
                ring2.bp, ring2.height))
            unvisited_edges.remove(edge)
            numxovers = math.gcd(ring2.numxovers, ring1.numxovers) - 1
            # calibrate and get first xover position
            refxoverpair = self.calibrate_ring(edge)
            # get all other xover positions
            xoverset = crossover.get_xovers(ring1, ring2, 'stap',
                                       threshold=edge.thresh, numxovers=numxovers)
            xoverset = crossover.distribute(xoverset, numxovers, 'stap')
            self.apply_crossovers(xoverset, 'stap')
        """
        For single wall structures
        """
        if not unvisited_edges.value:  # if list is empty
            log.out(__name__, "This seems to be a single wall routing.")
            log.out(__name__, "Visited all edges. Exiting staple crossover creation.")
            return  # exit early if single wall
        """
        For multiwall structures
        """
        log.out(__name__, "Staples along pathway complete. Connecting remaining edges.")
        
        for this_edge in unvisited_edges.getsorted():
            ring1 = this_edge.directed['from']
            ring2 = this_edge.directed['to']
            log.debug(__name__, "Connecting Edge( Ring({}), Ring({}) )"
                      .format(ring1.getsimpleposition(),
                              ring2.getsimpleposition()))

            numxovers = math.gcd(ring2.numxovers, ring1.numxovers) - 1

            # get all other xover positions
            xoverset = crossover.get_xovers(ring1, ring2, 'stap',
                                       threshold=edge.thresh, numxovers=numxovers)
            xoverset = crossover.distribute(xoverset, numxovers, 'stap')
            self.apply_crossovers(xoverset, 'stap')
        log.out(__name__, "Visited all edges. Exiting staple crossover creation.")

    def route_asym_staples(self):
        """
        Temporary function for testing methods specific to routing the asymmetric example 'Clover'
        Should be later consolidated into 'route_staples' method
        :return: None
        """
        log.out(__name__, "Proposed initial placement of staple crossovers.")
        # Initiate a crossover set
        self.xover_set = heuxover.CrossoverSet(self)
        # self.xover_set.init_basic(bondlen=config.INTERHELICAL + 0.03,  # That's about 1 nt offset left or right
        #                           xoverspacing1=config.VALIDXOVERSPACING_SAME,
        #                           xoverspacing2=config.VALIDXOVERSPACING_ADJ)


