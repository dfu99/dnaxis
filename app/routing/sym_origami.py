#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 14:52:27 2018

@author: dfu

CHANGELOG:
    backlog:
        - added connection pathway
        - added angle evaluation and connections via eulerian distance
        - added ring calibration to line up crossovers by scaffold reference
    7/17/2018:
        - Added scaffold connections
    7/27/2018:
        - backlog:
            - added staple crossovers
            - added scaffold breaking
        - added staple breaking
            - bugged, appears to get stuck in while loop of find_closest_and_break
            - pending, better strand length optimization algorithm
    8/12/2018:
        - bugfix
            - flipped true/false of isXover
    8/16/2018:
        - Did some housekeeping on section headers
        - Added some functions for identifying and working with regions
        - Complete new break_staples method by iterating through regions
    8/27/2018:
        - Added get_ringtoring_distance
        - Makes sense to have with get_ringtoring_angle
        - Helpful for debugging
        - Added debug_simple_adjacency
    8/31/2018:
        - Added multiwall compatibility for pathing
        - self.routing removed from more functions
    9/2/2018:
        - Added optional max distance argument to find_closest_and_break
    9/6/2018:
        - Bug fixed, endless looping in get_all_strand_regions
        - Initial condition crossover nucleotide sometimes on wrong side
        - Interfered with break condition
    9/12/2018:
        Adjacency threshold to eval_connections()
        Moved top level functions to their respective sub-files
        
"""

from .helper import log, mymath

# sub-classes
from . import (
        initsym,
        sequence,
        nicking,
        export,
        path,
        crossover,
        scxover,
        stxover,
        calibrate,
        marking,
        query,
        collide,
        scnicking,
        optimize,
        modules,
        validate,
        strand,
        heuxover,
        post,
        stats,
        functionalize)

from .score import (
        seeding)


class planes2grid:
    def __get__(self, instance, owner):
        return self.tolist(instance.planes)
    
    def __set__(self, instance, value):
        raise AttributeError("Value cannot be set.")
        
    def tolist(self, planes):
        grid = []
        for plane in planes:
            flat = []
            for module in plane.modules:
                flat.append(module)
            grid.append(flat)
        return grid


class Origami(sequence.Mixin,
              nicking.Mixin,
              export.Mixin,
              path.Mixin,
              marking.Mixin,
              crossover.Mixin,
              scxover.Mixin,
              stxover.Mixin,
              calibrate.Mixin,
              query.Mixin,
              collide.Mixin,
              scnicking.Mixin,
              seeding.Mixin,
              optimize.Mixin,
              validate.Mixin,
              strand.Mixin,
              heuxover.Mixin,
              post.Mixin,
              stats.Mixin,
              functionalize.Mixin):

    grid = planes2grid()

    def add_plane(self, plane):
        """
        Adds a plane to the Origami's plane list
        :param plane: see modules.Plane
        :return: None
        """
        self.planes.append(plane)
        plane.settop(self)

    def add_module(self, module):
        """
        Adds a module to the Origami and automatically assigns it to the correct plane
        :param module: see modules.Shape
        :return: None
        """
        # Check if we can add to existing Plane
        for p in self.planes:
            if p.height == module.height:
                p.add_module(module)
                # print("Adding to existing Plane found for Module height {}.".format(module.height))
                return
        # Otherwise add a new Plane
        newplane = modules.Plane(module.height)
        newplane.add_module(module)
        self.add_plane(newplane)
        # print("Adding new Plane {} for Module.".format(module.height))

    def __init__(self, rings, pathway, connections, filename, outputdir):
        """
        Initializes topological information
        :param rings: Module descriptors
        :param pathway: Routing of scaffold strand
        :param connections: Placement of all crossovers
        :param filename: Job name
        :param outputdir: Job location
        """
        self.filename = filename
        self.outputdir = outputdir
        self.man_path = pathway
        self.man_connections = connections
        self.extensions = None
        self.rings = rings
        motif = 's'
        log.out(__name__, "Creating ORIGAMI (filename =", filename, "; option =", motif, ")", headerlevel=0)

        # Keeps track of the highest numbered nucleotide index
        self.last_idx = 0  # Pending implementation
        # Keeps track of skipped crossover edges
        self.removed = path.EdgeSet([])

        # Keeps track of all merged strands
        self.strands = mymath.FakeSet([])

        # Keeps track of all Nicks
        self.nicks = []

        # Keeps track of GapNuclPairs
        # TODO: Get rid of this, this is horrible
        self.gapnuclpairs = []

    def generate(self):
        """
        The rendering process
        :return: None
        """
        self.load()
        self.route()
        self.save()

    def load(self):
        """
        Populate objects into the origami to create the template
        :return: None
        """
        self.planes = initsym.init_planes(self.rings)
        # Set the topology
        for plane in self.planes:
            plane.settop(self)

        # Initialize protected nucleotides set
        self.init_protected()

        # apply the connections file to link helices together
        self.apply_connections(self.man_connections)
        self.pathway = self.apply_pathway(self.man_path)
        self.pathway_edges = self.path_to_edges()
        self.pathway_strands = None

    def route(self):
        """
        Perform routing of crossovers and nicks
        :return: None
        """
        # Add staple crossovers
        self.route_staples()

        # Add scaffold crossovers
        self.route_scaffold()
        
        # Clean up colliding crossovers
        # possibly unnecessary to consolidate
        # self.clean()

        # Add nicks
        self.nick_all()

        # Apply sequence to scaffold
        self.apply_scaf_seq()
        self.optimize_seeds()
        log.out(__name__, "Structure generation complete.")

    def save(self):
        """
        Save outputs
        :return: None
        """
        log.out(__name__, "Saving outputs.")
        # Export to usable file types
        self.print_csv()  # staple sequences
        self.print_functionalized_csv('H', '53c')
        self.save_input()  # input rings
        self.save_conn()  # connections
        self.save_path()  # pathway
        self.print_oxDNA(print_crossovers=True)  # oxDNA conf and top
#        self.print_oxDNA_htrap(filename)
#         self.cdna_out = self.print_caDNAno(self.filename)  # cadnano json

    def renumrings(self):
        num = 1
        for r in self.get_modules():
            r.ringnum = num
            num += 1
