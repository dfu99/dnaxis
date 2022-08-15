#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scxover.py
Methods upon the scaffold crossovers
Created on Thu Oct 11 17:15:40 2018
As fork from crossover.py

Changelog:
    

@author: dfu

"""

from . import crossover, sequence, path, noncyclic
from app import config
from .helper import mymath, strandnav, log
import random


class Mixin:
    # =============================================================================
    # Scaffold stuff
    # =============================================================================

    # top level scaffold routing and crossovers function
    # iterate through connection pathway
    # rotate about routing 1/3 each time
    # connect closest scaffold crossover
    # TODO: There is maybe a wiser design criteria than 1/3 rotation
    def route_scaffold(self):
        log.out(__name__, "Adding scaffold crossovers to routing and removing colliding staple crossovers.")

        if self.is_multi_scaffold():
            log.out(__name__, "This origami seems to require multiple scaffolds. "
                              "Searching for best segm entation position.")
            removed = self.break_plane_asym()
            self.removed = removed
            for r in list(removed.value):
                log.debug(__name__, "Skip scaffold crossover {}.".format(r.toString()))
        else:
            log.out(__name__, "This origami seems to require only one scaffold. Segmentation skipped.")
            removed = path.EdgeSet([])
            self.removed = removed

        # set an initial position
        pos = 45  # TODO: Should depend on crossover motif 3, 4, or 5

        # iterate through pathway
        for edge in self.pathway_edges:
            if removed.contains(edge):
                log.out(__name__, "Structure segmented at Edge({}).".format(edge.toString()))
                continue
            ring1 = edge.directed['from']
            ring2 = edge.directed['to']
            log.debug(__name__, "Connecting scaffold for Ring({}, {}), Ring({}, {})."
                      .format(ring1.bp, ring1.height, ring2.bp, ring2.height))
            numxovers = 1  # Fixed for scaffold
            xovers = crossover.get_xovers(ring1, ring2, 'scaf',
                                          threshold=config.VALIDXOVERTHRESHBP + (edge.thresh - config.INTERHELICAL),
                                          numxovers=numxovers)
            # Get a distribution based on fixed crossover count
            xoverset = crossover.distribute(xovers, numxovers, 'scaf')
            self.apply_crossovers(xoverset, 'scaf')
            pos = (pos + 120) % 360  # TODO: Should depend on crossover motif 3, 4, or 5

    def is_multi_scaffold(self):
        """
        # return True if routing is larger than 1 scaffold
        # False otherwise
        :return:
        """
        log.out(__name__, "Checking if the structure requires multiple scaffolds.")
        size = self.get_scaf_size()
        numscafs = len(sequence.chooseseqs(size))
        if numscafs == 0:  # null case
            raise ValueError("ERROR: No scaffolds were applied.")
        elif numscafs == 1:
            return False
        else:
            return True

    def breakplane(self):
        """
        Return the edge along the pathway that should be skipped in a multiscaffold routing

        # Termination conditions:
        # 1. If scaffolds are used up before rings are traversed, raise error
        # 2. Once rings are fully traversed, return the removed list

        # Process:
        # 1. Start with the first scaffold
        # 2. Record its length

        # 3. Start with the first ring in the pathway
        # 4. Subtract the length of the ring from length of scaffold

        # 5. Move to the next ring and subtract again
        # 6. Repeat step 5
        # 6a. After each subtraction, check if the length of the scaffold is used
        # 6b. If it is, backtrack until we get out of the current, incompletely filled plane
        # 6c. create that edge and add it to the removed list

        # 7. Now try getting the next scaffold
        # 7a. Except StopIteration, raise an error that there is no scaffolds left
        # 7b. Otherwise get the length of the new scaffold
        # 8. Continue with traversing the rings as in step 3
        :return: Removed (list) of Edge of connections to skip in the pathway
        """

        # Get required scaffold length of entire origami structure
        size = self.get_scaf_size()
        # Get enough scaffolds to complete the routing
        scafs = sequence.chooseseqs(size)
        iter_scafs = iter(scafs)
        # Initialize custom set type EdgeSet for removed collection (output)
        removed = path.EdgeSet([])
        self.removed = removed
        # Convert the pathway graph of VERTICES denoting the scaffold routing pathway to a custom queue type Tape
        rings = mymath.Tape(self.pathway)

        while True:
            try:
                next_scaf = next(iter_scafs)
                avail = len(sequence.getseq(next_scaf))
                log.out(__name__, "Retrieved and applying {}, length {}.".format(next_scaf, avail))
            except StopIteration:
                raise RuntimeError("No scaffold remaining")
            while True:
                try:  # Get next ring
                    thisring = rings.forward()
                except StopIteration:  # Reached the last ring
                    return removed
                avail -= thisring.bp
                if avail < 0:  # Scaffold used up
                    log.out(__name__, "Used up {}. Getting next scaffold.".format(next_scaf))
                    while True:  # Step back to plane transition
                        thisring = rings.current()
                        lastring = rings.backward()

                        thisplane = thisring.up()
                        lastplane = lastring.up()
                        if thisplane != lastplane:  # different planes
                            # plan to remove this egde from scaffold xovers
                            # 20200220 EdgeSet was converted to List instead of Set
                            #   This is because set may add some randomness to the order
                            #   List does not have an add method
                            #   But this is preserved to debug later, as I do not want to go too far down this
                            #   rabbit hole for now.
                            removed.add(path.Edge(lastring, thisring))
                            break
                        else:  # still same plane
                            pass  # restarts loop at current ring
                    next_scaf = next(iter_scafs)
                    avail = len(sequence.getseq(next_scaf))
                    log.out(__name__, "Retrieved and applying {}, length {}.".format(next_scaf, avail))
                else:  # continue
                    pass

    def break_plane_asym(self):
        """
        Asymmetric version of break plane that uses edge(m1, me2) version of pathway instead of m1, m2, m3... sequence
        version of pathway
        :return: Removed (list) of Edge of connections to skip in the pathway
        """

        # Get required scaffold length of entire origami structure
        size = self.get_scaf_size()
        # Get enough scaffolds to complete the routing
        scafs = sequence.chooseseqs(size)
        iter_scafs = iter(scafs)
        # Initialize custom set type EdgeSet for removed collection (output)
        removed = path.EdgeSet([])
        # Convert the pathway graph of EDGES denoting the scaffold routing pathway to a custom queue type Tape
        edges = mymath.Tape(self.pathway)

        skip = [((64, 0), (65, 0)), ((64, 1), (65, 1)), ((63, 1), (66, 1))]
        # skip = []

        while True:
            try:  # Checks if there is another scaffold
                next_scaf = next(iter_scafs)
                avail = len(sequence.getseq(next_scaf))
                log.out(__name__, "Retrieved and applying {}, length {}.".format(next_scaf, avail))
            except StopIteration:
                raise RuntimeError("No scaffold remaining")

            # Run this loop if there are remaining scaffolds
            while True:
                try:  # Get next edge
                    thisedge = edges.forward()
                    fid = thisedge[0]  # Walk on the source vertex of each tuple
                    this_module = self.planes[fid[0]].modules[fid[1]]
                except StopIteration:  # Reached the last edge
                    return removed
                avail -= len(strandnav.absgetstrand(this_module.helix.scaf[0]))
                if not thisedge[1]:  # No target module for path, indicates it was manually terminated
                    try:
                        nextedge = edges.valueAt(edges._ind+1)
                        nid = nextedge[0]
                        next_module = self.planes[nid[0]].modules[nid[1]]
                        removed.add(path.Edge(this_module, next_module))
                        log.system("Detected a manual skip in pathway input between {} and {}. "
                                   "Adding it to the removed list ignoring remaining length of scaffold.".format(
                            fid, nid))
                        next_scaf = next(iter_scafs)
                        avail = len(sequence.getseq(next_scaf))
                        log.out(__name__, "Retrieved and applying {}, length {}.".format(next_scaf, avail))
                    except IndexError:  # Should have hit the end of the pathway
                        return removed
                elif thisedge in skip:
                    nid = thisedge[1]
                    next_module = self.planes[nid[0]].modules[nid[1]]
                    removed.add(path.Edge(this_module, next_module))
                    log.system("Hardcoded pathway skip between {} and {}. "
                               "Adding it to the removed list ignoring remaining length of scaffold.".format(fid, nid))
                    next_scaf = next(iter_scafs)
                    avail = len(sequence.getseq(next_scaf))
                    log.out(__name__, "Retrieved and applying {}, length {}.".format(next_scaf, avail))

                elif avail < 0:  # Scaffold used up
                    log.out(__name__, "Used up {}. Retrieving next scaffold.".format(next_scaf))
                    while True:  # Step back to plane transition
                        thisedge = edges.current()
                        lastedge = edges.backward()
                        fid = lastedge[0]
                        last_module = self.planes[fid[0]].modules[fid[1]]
                        tid = lastedge[1]
                        this_module = self.planes[tid[0]].modules[tid[1]]

                        # Option: modular_plane
                        # Ensures that every plane is fully filled with a single scaffold
                        if config.ROUTING == "modular_plane":
                            thisplane = this_module.up()
                            lastplane = last_module.up()
                            if thisplane != lastplane:  # different planes
                                # plan to remove this edge from scaffold xovers
                                removed.add(path.Edge(last_module, this_module))
                                break
                            else:  # still same plane
                                pass  # restarts loop at current ring

                        # Option: modular_strand
                        # Ensures that every helix is fully filled with a single scaffold (same plane can have multiple
                        # scaffolds
                        elif config.ROUTING == "modular_strand":
                            this_strand = this_module.helix.scaf[0].get_top_strand()
                            last_strand = last_module.helix.scaf[0].get_top_strand()
                            if this_strand != last_strand:  # Different helices
                                removed.add(path.Edge(last_module, this_module))
                                break
                            else:  # Still same strand (Shouldn't ever happen in "strand" case
                                pass  # Restarts loop at current ring
                        else:
                            raise NotImplementedError("There is no control logic to handle option {}."
                                                      .format(config.ROUTING))
                    next_scaf = next(iter_scafs)
                    avail = len(sequence.getseq(next_scaf))
                    log.out(__name__, "Retrieved and applying {}, length {}.".format(next_scaf, avail))
                else:  # continue
                    pass

    def route_asym_scaffold(self):
        """
        Temporary function for testing methods specific to routing the asymmetric example 'Clover'
        Should be later consolidated into 'route_scaffold' method
        :return: None
        """
        log.out(__name__, "Adding scaffold crossovers to routing and removing colliding staple crossovers.")

        '''Check configuration settings'''
        # If config/ROUTING was not set correctly.
        if config.ROUTING.lower() not in ['modular_plane', 'modular_strand']:
            raise RuntimeError("Invalid option {} for ROUTING. "
                               "Choose modular by plane or strand.".format(config.ROUTING))
        # Checks if the routing requires multiple scaffolds
        # Modular ensures a helix or plane is always fully filled with a single scaffold
        # If so, it gets a list of crossovers to skip
        elif self.is_multi_scaffold() and config.ROUTING.lower().startswith('modular'):
            log.out(__name__, "This origami seems to require multiple scaffolds. "
                              "Searching for best segmentation position.")
            removed = self.break_plane_asym()
            self.removed = removed
            for r in list(removed.value):
                log.out(__name__, "Skipping scaffold crossover {}.".format(r.toString()))
        # Single scaffold sets the removed list to an empty list
        else:
            log.out(__name__, "This origami seems to require only one scaffold."
                              "Segmentation skipped.")
            removed = path.EdgeSet([])
            self.removed = removed

        '''Find and add scaffold crossovers'''
        # Iterate through pathway and add single scaffold crossovers to link up all modules together
        # Save a flag in the rare case that a crossover was removed following a gap xover (which offsets the every
        #   other pattern that they should be placed by)
        gap_flag = False
        for edge in self.pathway_edges:
            module1 = edge.directed['from']
            module2 = edge.directed['to']
            # DEBUG
            # print("Scaffold connection on {} {} to {} {} using {}".format(module1.note, module1.get_top_plane().yaw,
            #                                                               module2.note, module2.get_top_plane().yaw,
            #                                                               edge.thresh))
            if removed.contains(edge):
                # Would have placed a gap half xover here
                if noncyclic.gap_exists(module1, module2) and not noncyclic.gap_xover_exists(module1, module2):
                    log.system("[WARNING]: An edge containing a gap scaffold crossover was skipped. "
                               "Leftover strand may be dangling.")
                    gap_flag = True
                log.out(__name__, "Structure segmented at Edge({}).".format(edge.toString()))
                continue

            if noncyclic.gap_exists(module1, module2) and \
                    (not noncyclic.gap_xover_exists(module1, module2) or gap_flag):
                # Two half crossovers placed at fix position determined by gap in structure
                xoverset = crossover.get_gap_xovers(module1, module2)
                gap_flag = False
            else:
                xoverset = crossover.get_xovers(module1, module2, 'scaf', threshold=config.VALIDXOVERTHRESHBP + (edge.thresh - config.INTERHELICAL))
                try:
                    xoverset = [random.choice(xoverset)]
                except IndexError:
                    raise IndexError("Cannot choose from empty sequence of crossovers between {} and {}".format(module1, module2))

            if not xoverset:
                raise RuntimeError("No crossovers were applied for edge between {} and {}.".format(module1, module2))
            # No need to rotate for asymmetric structures because each continuous segment is made of multiple segments
            # and crossover pathway routing can be consolidated to a certain region.
            self.apply_asym_scaffold_crossovers(xoverset)

    def apply_asym_scaffold_crossovers(self, xoverset):
        """
        Temporary function for testing methods specific to routing the asymmetric example 'Clover'
        :param xoverset: List of xoverpairs to apply
        :return: None
        """
        for xp in xoverset:
            fx = xp.apply()
            self.protect(xp.n1.to_tuple[0], config.SCPROTECT)
            self.protect(xp.n2.to_tuple[0], config.SCPROTECT)
            fx.settype('scaf')
            xp.n1.up().addxover(fx)
            xp.n2.up().addxover(fx)
