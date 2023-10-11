"""
Created by Daniel Fu (Reif Lab, Duke University) at 4/10/2021

Name        :
Project     :
Description :
Interpreter : Python 3.7.4
"""

import itertools
from app import config
from .exceptions import *
from .helper import mymath, strandnav, log, mytrig
import os
import numpy as np


class Mixin:
    # =============================================================================
    # stats
    # =============================================================================
    def num_xovers(self):
        """
        Counts the number of crossovers in the structure
        :return:
        """
        all_modules = self.get_modules()
        all_xovers = []
        for m in all_modules:
            for x in m.objxovers:
                all_xovers.append(x)
        while True:
            for x1, x2 in itertools.permutations(all_xovers, 2):
                if x1 == x2:
                    all_xovers.remove(x2)
                    break
            if all(not x1 == x2 for x1, x2 in itertools.permutations(all_xovers, 2)):
                break
        # for entry in all_xovers:
        #     log.debug(entry)
        return len(all_xovers)

    def num_scaf_xovers(self):
        count = 0
        for xover in self.get_all_xovers():
            if xover.strandtype == 'scaf':
                count += 1
        return count

    def num_stap_xovers(self):
        count = 0
        for xover in self.get_all_xovers():
            if xover.strandtype == 'stap':
                count += 1
        return count

    def num_xovers_by_jstrand(self, jstrand1, jstrand2):
        """
        Counts the number of crossovers between two strands
        :param js1:
        :param js2:
        :return:
        """
        ref_base1 = jstrand1.modules[0].helix.stap[0]
        ref_base2 = jstrand2.modules[0].helix.stap[0]

        bases_list1 = strandnav.absgetstrand(ref_base1)
        bases_list2 = strandnav.absgetstrand(ref_base2)
        # Initially begin counting at the latter half of a crossover
        # So that when coming back around, the distance between the first crossover and last crossover will
        # still be counted
        for idx in range(len(bases_list1)):
            this_base_idx = idx
            next_base_idx = (idx + 1) % len(bases_list1)
            this_base = bases_list1[this_base_idx]
            next_base = bases_list1[next_base_idx]
            # If a pair of nucleotides are both part of a crossover
            if ((strandnav.isxover(this_base) or strandnav.isxover(this_base.Comp)) and
                    (strandnav.isxover(next_base) or strandnav.isxover(next_base.Comp))):
                init_base_idx = next_base_idx
                init_base = bases_list1[init_base_idx]
                # Get either the scaffold or staple crossover
                try:
                    last_xover = self.findxover(init_base)
                except NuclNotFound:
                    last_xover = self.findxover(init_base.Comp)
                break
            else:
                init_base_idx = 0
                last_xover = None

        # Cycle the list of crossovers to begin at the latter half of it so it is not counted twice
        helix_xovers = mymath.cycle_list(bases_list1, init_base_idx)
        crossover_count = 0

        # Go through every base in the strand
        for base in helix_xovers:
            # Only check crossovers
            if strandnav.isxover(base) or strandnav.isxover(base.Comp):
                # checks for either staple or scaffold crossover
                try:
                    xobj = self.findxover(base)
                except NuclNotFound:
                    xobj = self.findxover(base.Comp)
                xover_nucls = [xobj.n1.n5, xobj.n1.n3, xobj.n2.n5, xobj.n2.n3]
                if xobj == last_xover:  # Ignore consecutive bases that are part of same crossover
                    pass
                elif (any(nucl in bases_list1 for nucl in xover_nucls) and
                      any(nucl in bases_list2 for nucl in xover_nucls)):  # Staple crossover
                    crossover_count += 1
                    last_xover = xobj  # Reset last reference point
                elif (any(nucl.Comp in bases_list1 for nucl in xover_nucls) and
                      any(nucl.Comp in bases_list2 for nucl in xover_nucls)):  # Scaffold crossover
                    crossover_count += 1
                    last_xover = xobj  # Reset last reference point
                else:
                    pass

        return crossover_count

    def print_long_bonds(self, dnaorigami, outputdir):
        f = open(os.path.join(outputdir, "bonds.txt"), 'w')
        all_xovers = dnaorigami.get_all_xovers()
        for xover in all_xovers:
            n1c = xover.n1.coords
            n2c = xover.n2.coords
            f.write("{}\nDistance:{}\n".format(xover, np.linalg.norm(n2c-n1c)))
        f.close()

    def crossover_spacing_by_jstrand(self, jstrand1, jstrand2):
        """
        Logs the spaces between crossovers spanning two strands
        :param jstrand1: JoinedStrand
        :param jstrand2: JoinedStrand
        :return:
        """
        spacings = []
        # Find an xover spanning the input joined strands to begin on
        for xover in jstrand1.modules[0].objxovers:
            if xover.spans(jstrand1, jstrand2) and xover.strandtype != 'scaf':
                init_xover = xover
                break
        # Get the NuclPair 5' end on jstrand1 as an initial position
        try:
            init_nucl = init_xover.onmodule(jstrand1.modules[0]).n5
        except UnboundLocalError:
            log.log_warning("Could not calculated crossover spacings per strand. "
                            "Something is probably wrong with the staple crossovers.")
            return [0]
        next_nucl = init_nucl.__strand5__
        dist = 2
        while True:
            # Determine if same edge crossover
            for xover in jstrand1.modules[0].objxovers:
                if next_nucl in xover and xover.spans(jstrand1, jstrand2):
                    next_xover = xover
                    next_nucl = next_xover.onmodule(jstrand1.modules[0]).n5
                    spacings.append(dist)
                    dist = 1
                    if next_xover == init_xover:
                        return spacings
                elif next_nucl.Comp in xover and xover.spans(jstrand1, jstrand2):
                    next_xover = xover.complement()
                    next_nucl = next_xover.onmodule(jstrand1.modules[0]).n5
                    spacings.append(dist)
                    dist = 1
                    if next_xover == init_xover:
                        return spacings
            next_nucl = next_nucl.__strand5__
            dist += 1

    def twist_spacing_by_jstrand(self, jstrand1, jstrand2):
        """
        Logs the number of base pairs traversed for an amount of turn
        Take the intended angle of the next strand as the true angle of any crossover
        :param jstrand1:
        :param jstrand2:
        :return:
        """
        twists = []
        # Find an xover spanning the input joined strands to begin on
        for xover in jstrand1.modules[0].objxovers:
            if xover.spans(jstrand1, jstrand2) and xover.strandtype != 'scaf':
                init_xover = xover
                break
        # Get the NuclPair 5' end on jstrand1 as an initial position
        init_nucl = init_xover.onmodule(jstrand1.modules[0]).n5
        this_theta = init_nucl.theta
        next_nucl = init_nucl.__strand5__
        next_theta = next_nucl.theta

        angle_sum = mytrig.acute_diff(this_theta, next_theta)
        # UNFINISHED
        pass
        # STUB

    def pap_xovers(self):
        """
        Shows whether a crossover is a parallel (False) or anti-parallel (True) crossover
        :return:
        """
        all_modules = self.get_modules()
        all_xovers = []
        for m in all_modules:
            for x in m.objxovers:
                all_xovers.append(x)
        while True:
            for x1, x2 in itertools.permutations(all_xovers, 2):
                if x1 == x2:
                    all_xovers.remove(x2)
                    break
            if all(not x1 == x2 for x1, x2 in itertools.permutations(all_xovers, 2)):
                break
        for entry in all_xovers:
            log.debug("Crossover orientation: {}".format(entry.is_anti_parallel()))

    def junction_xovers(self):
        """
        Collection of analytics for crossovers that span junctions between different scaffolds
        :return:
        """
        # Find each junction
        # Go through every connection and check if there is a scaffold crossover
        # If there isn't, save those edges in this list as scaffold junctions
        junctions = []
        for key in self.man_connections:
            module1 = self.planes[key[0]].modules[key[1]]
            for val in self.man_connections[key]:
                module2 = self.planes[val[0]].modules[val[1]]
                js1 = module1.get_top_strand()
                js2 = module2.get_top_strand()
                # If there is any FullXover of strandtype scaffold spanning the two JoinedStrands,
                #   then it is not a junction
                if not any(xobj.spans(js1, js2) and xobj.strandtype == 'scaf' for xobj in module1.objxovers):
                    # If any of the existing tuples have a permutation already in the list, don't add it
                    if not any(key in kv and val in kv for kv in junctions):
                        junctions.append((key, val))
        # print("Junctions are at", junctions)

        # For each junction, build another list of their crossovers.
        # Then for each of those crossovers, show the length of the strand.
        for kv in junctions:
            # print("At junction", kv)

            key = kv[0]
            val = kv[1]
            module1 = self.planes[key[0]].modules[key[1]]
            # print("k:", module1)
            module2 = self.planes[val[0]].modules[val[1]]
            # print("v:", module2)
            js1 = module1.get_top_strand()
            js2 = module2.get_top_strand()
            xoverlist = []
            # Query the xovers at the junction
            for xobj in module1.objxovers:
                if xobj.spans(js1, js2):
                    xoverlist.append(xobj)

            # # Try and automatically extend the staple to adjacent strands
            # for xobj in xoverlist:
            #     strand1 = strandnav.getstrand(xobj.n1.n5)
            #
            #     print("Strand len: {} helices: {} order: {}".format(len(strand1), count_jumps(strand1)+1,
            #                                                         self.strand_helices(strand1)))
            #     # Try and automatically extend the staple to adjacent strands
            #     # end5 = strand1[0]
            #     # end3 = strand1[-1]
            #     # adj5 = strandnav.getstrand(end5.__strand5__)
            #     # adj3 = strandnav.getstrand(end3.__strand3__)
            #     # print("If strand is extended")
            #     strand2 = strandnav.getstrand(xobj.n1.n3)
            #     print("Strand len: {} helices: {} order: {}".format(len(strand2), count_jumps(strand2) + 1,
            #                                                         self.strand_helices(strand2)))

    def strand_helices(self, strand):
        """
        Returns the strand as a sequential list of helices it traverses
        :param strand: list of Nucleotide
        :return: list of module index tuples
        """
        order = []
        for n in strand:
            m = n.get_top_module()
            ind = self.get_module_index(m)
            if ind not in order:
                order.append(ind)
        return order


def count_jumps(strand):
    """
    Counts the number of helices the strand traverses
    :param strand: list of Nucleotide
    :return: int
    """
    jumps = 0
    for n1, n2 in zip(strand[:-1], strand[1:]):
        if n1.get_top_strand() != n2.get_top_strand():
            jumps += 1
    return jumps


if __name__ == "__main__":
    pass
