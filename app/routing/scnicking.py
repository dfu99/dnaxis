#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  1 21:24:23 2019

@author: dfu

scnicking.py
"""
from .helper import strandnav, log, mymath
from . import nicking, crossover
from app import config


class Mixin:
    # break the scaffold strand
    # start anywhere and look for a clear regions
    def nick_scaffold(self):
        """
        This is the automatic scaffold nicking algorithm
        It proceeds by the following algorithm:
            1. Finds all the partitions
            2. Then matches these partitions to sections of the scaffold routing pathway
            3. For the first and last partitions of the Origami, places the nick on the first or last helix (facing the
            outside of the structure)
            4. For middle partitions, it doesn't really matter, but place the nick on the first helix of the section
        :return: List of placed nicks for recordkeeping
        """
        # Get the partitioned sections of the Origami from where their scaffold crossovers have been skipped
        partitions = self.partition_from_removed(self.removed)
        numnicks = len(partitions)
        log.out(__name__, "Looking for {} scaffold strand nick locations.".format(numnicks))
        scnicks = []  # List of nicks that will be placed
        # Consolidate the sequence of the pathway that is present in the partitioned section
        pathway_nodes = mymath.edgelist2nodelist(self.pathway_strands)
        for secnum, section_modules in enumerate(partitions):
            # Convert the modules to the JoinedStrands
            section_jstrands = [m.helix.scaf[0].get_top_strand() for m in section_modules]
            # Match the nodes in the pathway to nodes in the sections
            section_pathway_nodes = []
            for node in pathway_nodes:
                if node in section_jstrands:
                    section_pathway_nodes.append(node)

            # Start on the first node of the pathway if first section
            # Start on the last node of the pathway if last section
            # Start on the first node of the pathway if middle section
            if secnum != len(partitions):
                init_nucl = section_modules[0].helix.scaf[0]
            else:
                init_nucl = section_modules[-1].helix.scaf[0]
            nucl = init_nucl
            cycles = 0
            # Only valid if scaffold is not nicked yet
            if strandnav.checkloop(nucl):
                while True:
                    # Case 1, first pass. look for unprotected region
                    if cycles == 0:
                        if ((nucl in self.protected and nucl.__strand3__ not in self.protected) or
                                (nucl not in self.protected and nucl.__strand3__ in self.protected)):
                            break
                    # Case 2, subsequent passes.
                    if cycles > 0:
                        if strandnav.distfromfeature(nucl) <= config.PROTECT - cycles:
                            break
                    nucl = nucl.__strand3__
                    # Case 3, fully looped
                    if nucl == init_nucl:
                        log.debug(__name__, "Warning: No valid scaffold nick position" +
                                  "found for PROTECT={}. ".format(config.PROTECT - cycles) +
                                  "Reducing acknowledged protection range " +
                                  "to {}.".format(config.PROTECT - cycles - 1))
                        cycles += 1
                    # case 4, still haven't found anything good
                    if cycles >= config.PROTECT - config.LIMSPACINGOFFSET:
                        raise RuntimeError("Error: Was not able to place a scaffold nick anywhere reasonable.")

                self.protect(nucl, config.PROTECT)
                scnicks.append(nicking.Nick(crossover.NuclPair(nucl, nucl.__strand3__)))
            else:
                log.out(__name__, "WARNING: Tried to break scaffold loop but scaffold is already broken.")
        return scnicks  # returns the nicks that were created

    def scaffold_nick_on_ring(self, module):
        """
        Tries to place a nick on the input ring only
        :param module: Module(object)
        :return: List of Nick(object)
        """
        log.out(__name__, "Attempting to place a scaffold nick on module {}.".format(module))
        scnicks = []
        # Generaliation between symmetric and asymmetric structure
        # For an asymmetric structure, *.helix.scaf may not be the entire strand
        full_strand_scaf = strandnav.absgetstrand(module.helix.scaf[0])
        log.out(__name__, "Trying distance search method.")
        best_location = {"nucl": None, "dist": None}
        for nucl in full_strand_scaf:
            dist = strandnav.distfromfeature(nucl)
            if not best_location['nucl']:
                best_location['nucl'] = nucl
                best_location['dist'] = dist
            elif dist < best_location["dist"]:
                best_location['nucl'] = nucl
                best_location['dist'] = dist

        nucl = best_location['nucl']
        self.protect(nucl, config.PROTECT)
        scnicks.append(nicking.Nick(crossover.NuclPair(nucl, nucl.__strand3__)))
        return scnicks
