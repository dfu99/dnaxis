#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  1 21:24:23 2019

@author: dfu

scnicking.py
"""
from .helper import strandnav, log, mymath
from . import nicking, motifs
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
        for secnum, section_modules in enumerate(partitions):
            init_position_offset = 0
            while True:
                try:
                    # First, try position near the edges of the modules
                    if secnum != len(partitions):
                        init_nucl = section_modules[0+init_position_offset].helix.scaf[0]
                    else:
                        init_nucl = section_modules[-1-init_position_offset].helix.scaf[0]
                    nucl = init_nucl
                    # Only valid if scaffold is not nicked yet

                    assert strandnav.checkloop(nucl), log.out(__name__, "Tried to nick the scaffold but it already nicked. "
                                                                        "Something may have not run in the correct order.")

                    while True:
                        # Ideally we find a location that is not already marked as too close to another
                        # crossover or nick.
                        if any(n not in self.protected for n in (nucl, nucl.__strand3__)):
                            break
                        else:
                            # Go to next nucleotide
                            nucl = nucl.__strand3__
                        # If we've gone all the way around, fail, and go to the next module
                        if nucl == init_nucl:
                            raise RuntimeError
                    self.protect(nucl, config.PROTECT)
                    scnicks.append(nicking.Nick(motifs.NuclPair(nucl, nucl.__strand3__)))
                    break
                except RuntimeError:
                    init_position_offset += 1
                except IndexError:
                    raise SystemError("Was not able to place a nick anywhere reasonable on the structure.")
        return scnicks  # returns the nicks that were created

    def scaffold_nick_on_ring(self, module):
        """
        Tries to place a nick on the input ring only
        :param module: Module(object)
        :return: List of Nick(object)s
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
        scnicks.append(nicking.Nick(motifs.NuclPair(nucl, nucl.__strand3__)))
        return scnicks
