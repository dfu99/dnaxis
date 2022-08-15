#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created by Daniel Fu (Reif Lab, Duke University) on Sun Sep  9 19:20:19 2018

Name        : routing/query.py
Project     : cadaxisdna
Description : Queries that look through the entire DNA origami object
Interpreter : Python 3.7.4

Changelog:
    9/9/2018:
        File created
    9/27/2018:
        Added get_edges method
"""
import itertools
from . import path
from .helper import strandnav, mymath, log
from .exceptions import NuclNotFound


def remove_duplicate_strands(strandlist):
    """
    Goes through a list of strands and removes duplicates
    :param strandlist: List of list of nucleotides
    :return:
    """
    while True:
        for s1, s2 in itertools.permutations(strandlist, 2):
            if strandnav.strandequals(s1, s2):
                strandlist.remove(s2)
                log.log_warning("Duplicated strands were detected and removed.")
                break
        if all(not strandnav.strandequals(s1, s2) for s1, s2 in itertools.permutations(strandlist, 2)):
            break
    return strandlist


def get_xovers_btwn_rings(ring1, ring2):
    """
    Queries a list of crossovers that span ring1 and ring2
    :param ring1: address of ring1
    :param ring2: address of ring2
    :return: list of FullXover
    """
    xoverslist = {ring1: ring1.appliedxoversto(ring2), ring2: ring2.appliedxoversto(ring1)}
    return xoverslist


class Mixin:
    def get_all_regions(self):
        """
        Return all staple strands as regions
        :return:
        """
        all_staples = self.get_all_staples()
        all_regions = []
        next_region = []
        for staple in all_staples:
            for nucl in staple:
                if not next_region:
                    next_region.append(nucl)
                else:
                    if nucl.toFive != nucl.__strand5__:
                        all_regions.append(next_region)
                        next_region = [nucl]
        return all_regions


    def get_modules(self):
        """
        Returns all the modules in the structure in a 1D list
        20180825
            In singlewall structures, it was safe to loop once through self.routing.
            However in multiwall, we must look within each plane to get all the rings.
            It will be easier to keep the existing routing of eval_connections if all the rings are produced and we
            still iterate through a 1D list.
        :return: 1D list of all modules as objects
        """
        all_modules = []
        for plane in self.planes:
            for module in plane.modules:
                all_modules.append(module)
        return all_modules

    def get_module_index(self, module):
        """
        Returns the (plane, module) index associated with a module
        :param module: the module object
        :return: tuple index
        """
        for pid in range(len(self.planes)):
            for mid in range(len(self.planes[pid].modules)):
                if module is self.planes[pid].modules[mid]:
                    return pid, mid
        raise IndexError("Didn't find that module.")

    def get_module_by_index(self, ind):
        """
        Returns the module associated with an index (plane, module)
        :param ind: 2-tuple (plane index, module index)
        :return: DNA helix shape module object
        """
        try:
            pidx = ind[0]
            midx = ind[1]
        except IndexError as e:
            log.log_error("{} not valid module query. Use a (plane index, module index) 2-tuple".format(e, ind))

        try:
            return self.planes[pidx].modules[midx]
        except IndexError as e:
            log.log_error("Did not find module {}.".format(ind))

    def get_nick(self, nucl):
        """
        Returns the nick corresponding to input nucleotide ID
        :return: Nick
        """
        for n in self.nicks:
            if nucl.numid in n.nuclpair.numid:
                return n
        # print("Could not find a nick for this nucleotide. Context:")
        # print(nucl.numid)
        # print(nucl.toThree)
        # print(nucl.toFive)
        # print(nucl.__strand3__.toFive)
        # print(nucl.__strand5__)
        raise ValueError("Could not find nick for {}.".format(nucl.numid))

    def get_scaf_size(self):
        """
        Get scaf size adds up all the bp circumferences of all rings and returns the
        size of the origami routing
        """
        size = 0
        for plane in self.planes:
            for module in plane.modules:
                size += module.bp
        return size
    
    # get scaf size returns the size of a single helix
    def get_size(self):
        return self.get_scaf_size()*2
    
    # 20180927
    # Converts the adjacency per each module to a set of edges
    def get_edges(self):
        all_edges = []
        all_modules = self.get_modules()
        for module in all_modules:
            for adj_module in module.adjacent:
                all_edges.append(path.Edge(module, adj_module))
        return path.EdgeSet(all_edges)

    # returns both staples and scaffolds
    def get_all_strands(self):
        all_strands = self.get_all_staples() + self.get_all_scaffolds()
        all_strands = remove_duplicate_strands(all_strands)
        return all_strands
    
    # iterates through the origami and gets all the staple strands
    def get_all_staples(self):
        visited = []
        strands = []
        for plane in self.planes:
            for module in plane.modules:
                for nucl in module.helix.stap:
                    if nucl not in visited:
                        s = strandnav.getstrand(nucl)
                        if s not in strands:
                            strands.append(s)
                        for n in s:
                            visited.append(n)
        return strands
    
    # iterates through only the selected rings and gets the staple strands there
    def get_ring_staples(self, selection):
        visited = []
        strands = []
        for plane in self.planes:
            for module in plane.modules:
                if module in selection:
                    for nucl in module.helix.stap:
                        if nucl not in visited:
                            s = strandnav.getstrand(nucl)
                            strands.append(s)
                            for n in s:
                                visited.append(n)
        return strands
    
    # iterates through the origami and gets all the scaffold strands
    def get_all_scaffolds(self):
        visited = []
        strands = []
        for plane in self.planes:
            for module in plane.modules:
                for nucl in module.helix.scaf:
                    if nucl not in visited:
                        s = strandnav.getstrand(nucl)
                        if s not in strands:
                            strands.append(s)
                        for n in s:
                            visited.append(n)
        return strands

    def get_all_strand_loops(self, strandtype='stap'):
        """
        Go through all modules of the origami and consolidate loop motifs. Loop motifs are strands that are cyclic.
        A cyclic strand is one where there are 5' or 3' connections between the first [0] and last [-1] elements.

        If only this query is run, assumes that all strands in the object are loops
            or, that there are no nicks in the structure
        :return: List of loop strands
        """
        loops = []  # List of loops
        strands = []  # List of all other strands
        visited = []  # Visited nucleotides to avoid double counting
        all_modules = self.get_modules()
        # Visit all staple nucloetides by module
        for module in all_modules:
            if strandtype == 'scaf':
                tar_strand = module.helix.scaf
            else:
                tar_strand = module.helix.stap
            for nucl in tar_strand:
                if nucl in visited:  # Skip nucleotides already in a strand
                    pass
                else:
                    if strandnav.checkloop(nucl):  # TODO: This isn't very efficient
                        # 20181119 Saved as tuple rather than set to avoid randomness
                        loop = tuple(strandnav.getstrand(nucl))
                        # Add nucleotides of loop as visited
                        for l_nucl in loop:
                            visited.append(l_nucl)
                        loops.append(loop)
                    else:
                        strand = tuple(strandnav.getstrand(nucl))
                        for s_nucl in strand:
                            visited.append(s_nucl)
                        strands.append(strand)
                        # log.out(__name__, "Warning: There are nicks present in the structure already. Both loops and "
                        #                   "strands will be processed in the nicking algorithm. This may scramble "
                        #                   "intentional nicks.")
        return loops, strands

    def get_nucl(self, numid):
        """
        Returns the Nucleotide at the numid
        Directly queries the helix nucleotides
        :param numid: type int; id value of nucleotide
        :return: Nucleotide object
        """
        all_nucl = self.get_all_helix_nucl()
        for nucl in all_nucl:
            if nucl.numid == numid:
                return nucl
        raise ValueError("Did not find that nucleotide.")

    def get_all_helix_nucl(self):
        """
        Gets all nucleotides in the structure by querying the helices.
        Tries to avoid any issues with strand or ring topology.
        :return: list of all Nucleotide(object)
        """
        all_nucl = []
        for p in self.planes:
            for m in p.modules:
                for nucl in m.helix.scaf:
                    all_nucl.append(nucl)
                for nucl in m.helix.stap:
                    all_nucl.append(nucl)
        return all_nucl

    def get_all_strand_nucl(self):
        """
        Gets all the nucleotides in the structure by querying the strands
        :return: list of all Nucleotide(object)
        """
        strands = self.get_all_strands()
        strand_nucl = []
        for s in strands:
            for nucl in s:
                strand_nucl.append(nucl)
        return strand_nucl

    def current_index(self):
        """
        Returns the highest used nucleotide index
        :return: Index attribute of origami
        """
        return self.last_idx

    def findxover(self, nucl):
        """
        Find the crossover object associated with nucl
        :param nucl: Nucleotide object
        :return: FullXover object
        """
        all_modules = self.get_modules()
        for module in all_modules:
            for xobj in module.objxovers:
                if any([xobj.n1.n5.numid == nucl.numid, xobj.n1.n3.numid == nucl.numid,
                        xobj.n2.n5.numid == nucl.numid, xobj.n2.n3.numid == nucl.numid]):
                    return xobj
        raise NuclNotFound

    def get_all_xovers(self):
        """
        Returns all the FullXover objects in the structure
        :return:
        """
        all_modules = self.get_modules()
        all_xovers = mymath.FakeSet([])
        for module in all_modules:
            [all_xovers.add(xobj) for xobj in module.objxovers]
        return all_xovers

