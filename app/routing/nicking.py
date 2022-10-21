#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  8 17:28:13 2018

@author: dfu

nicking.py

CHANGELOG:
    - 9/8/2018:
        - File created. Separated nicking functions from origami.py
    - 9/12/2018:
        Resolved runtime error in find_closest_and_break
        Added top level functions from sym_origami.py
    11/18/2018:
        version 2 instead of the really basic greedy algorithm from before
"""

from . import crossover, sequence
from .exceptions import *
from app import config

from .helper import log, dnaconnector, mymath, strandnav


class Nick:
    """
    Nick object remembers the NuclPair location it is placed at
    It is necessary to use this object such that it can be undone

    The alternative is to get strands, then check the 5' and 3' absolute adjacent nucleotides
    """
    def __init__(self, np, active=True):
        self.nuclpair = np
        if active:
            self.nick()

    def nick(self):
        log.debug(__name__, "Created nick on {}".format(self.nuclpair.numid))
        dnaconnector.nuclbreak(self.nuclpair.to_tuple[0])

    def undo(self):
        (n5, n3) = self.nuclpair.to_tuple
        log.debug(__name__, "Undo nick on {}".format(self.nuclpair.numid))
        dnaconnector.nuclconnect(n5, n3)

    def __repr__(self):
        return "Nick[{}, {}]".format(self.nuclpair.n5, self.nuclpair.n3)


class NickSet:
    """
    A collection of all the Nicks for an instance of the Origami
    """
    def __init__(self, dnaorigami):
        pass


class Mixin:
    """
    Top level break function
       Nicks all strands
    """
    def nick_all(self):
        log.out(__name__, "Running nicking algorithm.")
        # Automatic staple nicking
        self.nick_staples()
        # Automatic scaffold nicking
        if config.SCAF_NICKING == 'auto':
            sc_nicks = self.nick_scaffold()
            self.sc_starts = [n.nuclpair.n3 for n in sc_nicks]

        # Hard coded manually selected scaffold nick locations
        # Symmetric only for 2 scaffolds
        elif config.SCAF_NICKING == 'symmetric_double':
            firstnick = self.pathway[0]
            lastnick = self.pathway[-1]
            self.scaffold_nick_on_ring(firstnick)
            self.scaffold_nick_on_ring(lastnick)

        # Asymmetric uses edges
        elif config.SCAF_NICKING == 'asymmetric_single':
            ###
            # Manual method
            idx1 = self.pathway[-1][1]
            firstnick = self.planes[idx1[0]].modules[idx1[1]]
            self.scaffold_nick_on_ring(firstnick)

            ###
            # Tries to find the highest distance from any other feature
            # for n in self.get_all_helix_nucl():
            #     ndist = strandnav.distfromfeature(n)
            #     try:
            #         if strandnav.distfromfeature(n) < ndistmin:
            #             ndistmin = ndist
            #             nmin = n
            #     except NameError:
            #         ndistmin = ndist
            #         nmin = n
            # _ = nicking.Nick(crossover.NuclPair(nmin, nmin.__strand3__))

        elif config.SCAF_NICKING == 'asymmetric_double':
            idx1 = self.pathway[0][0]
            idx2 = self.pathway[-1][1]
            firstnick = self.planes[idx1[0]].modules[idx1[1]]
            lastnick = self.planes[idx2[0]].modules[idx2[1]]
            self.scaffold_nick_on_ring(firstnick)
            self.scaffold_nick_on_ring(lastnick)
        else:
            log.log_error("{} is not a valid option for {}.".format(config.SCAF_NICKING, 'SCAF_NICKING'))

        # Contingency nicking staple loops
        if config.LOOP_CHECK:
            stap_loops, _ = self.get_all_strand_loops(strandtype='stap')
            for loop in stap_loops:
                refnucl = loop[0]
                _ = break_loop(refnucl)
    """
    Nicking sub-functions
    """
    def nick_staples(self):
        """
        Staple nicking top level function
        :return: None
        """
        log.dev(__name__, "Getting all loops.")
        (loops, strands) = self.get_all_strand_loops()
        # Get the loops and strands. For any cyclic structure, there should be no nicks or strands yet. "strands" should
        # be empty. But if intentional nicks were placed before this step, then strands created by those nicks will also
        # be included in the algorithm, so nicks may get scrambled.
        log.dev(__name__, "Got {} loops.".format(len(loops)))
        nickmap = {}  # key: loop (list of nucl); value: Nicks on loop
        log.out(__name__, "Nicking staple into short strands.")
        nickmap = self.chop_loops(loops, nickmap)
        nickmap = self.chop_strands(strands, nickmap)
        nickmap = self.merge_on_nicks(nickmap)
        # TODO: Get rid of this, this is horrible, has to do with manually breaking and reconnecting for acyclic
        # print("[DEBUG nickmap NoneType]:")
        # for key in nickmap:
        #     print(key, nickmap[key])
        self.nicks = [n for sub in nickmap.values() for n in sub]

    def chop_loops(self, loops, nickmap):
        """
        Add nicks into looped strand to chop it up into a lot of small strands with length determined by NICKSPACING
        :param loops: List of loop strands
        :param nickmap: Dict mapping Nick to strand
        :return: Updated nickmap
        """
        for loop in loops:
            # Debugging
            # log.dev(__name__, "Loop:", headerlevel=4)
            # self.debug_strand(loop)
            '''
            # Consolidate substrands as protected or unprotected
            # @ invariant:
            #   even indices: protected
            #   odd indices: unprotected
            '''
            c_loop = self.split2substrands(loop)
            # Troubleshooting
            # for c in c_loop:
            #     log.dev(__name__,"Consolidated substrands in loop:")
            #     self.debug_strand(c)
            # In the following, we "chop up" the strand by creating many nicks
            if len(c_loop) == 1:
                # Only 1 substrand indicates a fully protected loop that cannot be excised
                nickmap[loop] = self.heuristic_nick(c_loop[0])
            else:
                # Add nicks to the ends of each substrand
                # log.dev(__name__, "Excising substrands.")
                nickmap[loop] = excise(c_loop)
                # log.dev(__name__, "Substrands excised.")
                # Then go through each substrand and add nicks
                for i, sub in enumerate(c_loop, start=0):
                    # If odd index, (unprotected)
                    if i % 2:  # i % 2 == 1
                        # log.dev(__name__, "Running normal nicking on strand {} of {}, with NICKSPACING={}".format(
                        #     i + 1, len(c_loop), config.NICKSPACING))
                        # self.debug_strand(sub)
                        nickmap[loop] += self.nickstrand(sub, config.NICKSPACING, ignoremarking=False)
                    # If even index (protected substrand)
                    else:  # i % 2 == 0
                        # log.dev(__name__, "Running forced nicking on strand {} of {}".format(i + 1, len(c_loop)))
                        # self.debug_strand(sub)
                        # if too long, iterate down from NICKSPACING until lengths are acceptable
                        if len(sub) > config.LENUP:
                            # Substrand is long and fully protected. Place nicks that may violate the spacing
                            # log.dev(__name__, "Substrand is long ({}) but protected. Calling heuristic_nick.".format(len(sub)))
                            # for n in sub:
                            #     log.dev(__name__, n.up().up().bp, n.up().up().height)
                            nickmap[loop] += self.heuristic_nick(sub)
                        # Otherwise leave it as is
                        else:
                            # Substrand is protected, but of acceptable length so take no action.
                            # log.dev(__name__, "Substrand is protected, but of acceptable length={}. "
                            #                   "No action taken.".format(len(sub)))
                            pass
        return nickmap

    def chop_strands(self, strands, nickmap):
        """
        Add nicks into strand to chop it up into a lot of small strands with length determined by NICKSPACING
        :param strands: List of strands to chop up
        :param nickmap: Dict mapping Nick to strand
        :return: Updated nickmap
        """
        for strand in strands:
            nickmap[strand] = self.nickstrand(strand, config.NICKSPACING, ignoremarking=False)
        return nickmap

    def merge_on_nicks(self, nick_map):
        """
        After having all the nicks, run a merging algorithm to merge short strands with adjacent strands, as long as a
        strand exceeding the upper limit is not formed.
        :param nick_map: dict linking strands to list of Nicks
        :return: None
        """
        log.out(__name__, "Running merger algorithm to combine short strands.")
        for strand_key in nick_map:
            # log.dev(__name__, "{} nicks for strand_key are:".format(len(nick_map[strand_key])))
            # for _n in nick_map[strand_key]:
            #     log.dev(__name__, "obj Nick: {}".format(_n.nuclpair.numid))
            # Remove unnecessary nicks through merge algorithm and update nickmap
            nicks = nick_map[strand_key] = self.merge_short(nick_map[strand_key])
            # log.dev(__name__, "Loop merge completed.")
            # Create new protected regions for those nicks
            self.protect_nicks(nicks, config.PROTECT)
            # log.dev(__name__, "Loop nick protection completed.")
        return nick_map

    def heuristic_nick(self, strand):
        """
        Heuristic process that places nicks according to a decrementing NICKSPACING until an acceptable set of
            substrands is created
        :param strand:
        :return:
        """
        # log.dev(__name__, "Started a heuristic nick on strand:")
        # self.debug_strand(strand)
        # upper bound NICKSPACING
        # lower bound 3
        for i in range(config.NICKSPACING, 2, -1):
            # log.dev(__name__, "Attempting NICKSPACING={}".format(i))
            _nicks = self.nickstrand(strand, i, ignoremarking=True)

            # Skip if nickstrand didn't add any nicks
            if _nicks:
                _strands = nickstostrands(_nicks)
                # Debugging
                # log.dev(__name__, "Produced {} strands from nicks".format(len(_strands)))
                # for _s in _strands:
                #     log.dev(__name__, "Strand:")
                #     self.debug_strand(_s)
                if any(len(_s) > config.LENUP for _s in _strands):  # If any strand is still too long, retry the loop
                    # log.dev(__name__, "Strand len>{} detected. Redoing...".format(config.LENUP))
                    for _n in _nicks:
                        # log.dev(__name__, "Undoing {} from {}".format(_n, _nicks))
                        _n.undo()
                else:  # All strands are of valid length
                    # lens = [len(s) for s in _strands]
                    # log.dev(__name__,
                    #         "Current strand slicing accepted. Shortest:{} bps; Longest:{} bps".format(min(lens),
                    #                                                                                   max(lens)))
                    return _nicks
            else:
                # The strand wasn't nicked at all, so we can skip the other calculations for this loop
                pass
                # log.dev(__name__, "No nicks aded for NICKSPACING={}".format(i))
        # Since heuristic_nick is the last resort nicking strategy, if it fails, the process should not continue
        raise HnickError("heuristic_nick error: no valid NICKSPACING was found.")

    def protect_nicks(self, nicks, num):
        if not nicks:
            return
        for nick in nicks:
            # start with the (5',3') pair
            nucl3 = nick.nuclpair.n3
            nucl5 = nick.nuclpair.n5
            self.add_protected(nucl3)
            self.add_protected(nucl5)
            self.add_protected(nucl3.Comp)
            self.add_protected(nucl5.Comp)
            # move outward num-1 spaces, such that the total number of nucleotides
            # visited from each 5' or 3' nucleotide outward is equal to num
            for i in range(num - 1):
                if nucl3.toThree == -1:
                    pass
                else:
                    nucl3 = nucl3.toThree
                if nucl5.toFive == -1:
                    pass
                else:
                    nucl5 = nucl5.toFive
                self.add_protected(nucl3)
                self.add_protected(nucl5)
                self.add_protected(nucl3.Comp)
                self.add_protected(nucl5.Comp)

    def merge_short(self, nicks):
        """
        Merge short strands to shortest adjacent strand until strands are within the lower and upper length bounds
        :param nicks: Nicks applied into the origami
        :return: Updated list of Nicks
        """
        if not nicks:  # Nicks is empty, nothing to do
            return nicks
        kill_counter = 0  # Track while loop steps
        invalid = mymath.FakeSet([])  # strands that i've tried to join but didn't have any valid adjacent strands
        while True:
            '''CHECK Loop termination conditions'''
            kill_counter += 1
            if kill_counter > 50000:
                # Loop termination condition
                #   A very large number to terminate the WHILE if it gets stuck
                raise RuntimeError("Merge algorithm hasn't terminated after a very long time. Killing the process.")
            strands = nickstostrands(nicks) - invalid  # Updated remaining strands
            if len(strands) == 0:
                # Loop terimnation conditions
                #   If no more strands can be merged
                break
            smin = tuple(mymath.lstlenmin(strands))  # Get the shortest strand
            # if len(smin) > config.LENLOW:
            #     # Loop termination condition
            #     #   If the shortest strand is also above the lower bound length
            #     break

            '''Perform merge if not terminated'''
            # 5', 3' ends of strand
            sn5 = strandnav.goto5(smin[0])
            sn3 = strandnav.goto3(smin[0])

            '''Get the strands 5', 3' adjacent strands'''
            # If it is not a cyclic strand, the endpoints will not have an adjacent nucleotide
            # 5', 3' adjacent strands
            # POSSIBLY DEPRECATED IF WE MAKE ACYCLIC STRANDS STILL HAVE AN ABSOLUTE CONNECTION TO THEIR OTHER ENDPOINT
            try:
                adj5 = strandnav.getstrand(sn5.__strand5__)
            except AttributeError:
                adj5 = None
            try:
                adj3 = strandnav.getstrand(sn3.__strand3__)
            except AttributeError:
                adj3 = None

            '''Check each case determine whether a nick can be merged'''
            # Gradually remove strands (add to invalid set) if their case is terminal
            # Case: Strand forms a loop
            if strandnav.seqequals(adj5, smin) or strandnav.seqequals(adj3, smin):
                invalid.add(smin)

            # Case: Handles above non-cyclic errors for the remaining 3' or 5' strand
            elif adj5 is None:  # Handles the AttributeError above, 3' is default shorter
                if len(adj3) + len(smin) > config.LENUP:  # Strand would be longer than the upper bound
                    invalid.add(smin)
                else:  # Valid
                    nicks = _merge_strands(sn3, sn3.__strand3__, nicks)
            elif adj3 is None:  # Handles the AttributeError above, 5' is default shorter
                if len(adj5) + len(smin) > config.LENUP:  # Strand would be longer than the upper bound
                    invalid.add(smin)
                else:  # Valid length
                    nicks = _merge_strands(sn5.__strand5__, sn5, nicks)

            # Case: General case
            # Greedy implementation to eliminate shortest strands, pick the shorter adjacent strand
            elif len(adj5) <= len(adj3):  # 5' adjacent strand is shorter
                if len(adj5) + len(smin) > config.LENUP:  # Strand would be longer than the upper bound
                    invalid.add(smin)
                else:  # Valid length
                    nicks = _merge_strands(sn5.__strand5__, sn5, nicks)
            elif len(adj5) > len(adj3):  # 3' adjacent strand is shorter
                if len(adj3) + len(smin) > config.LENUP:  # Strand would be longer than the upper bound
                    invalid.add(smin)
                else:  # Valid
                    nicks = _merge_strands(sn3, sn3.__strand3__, nicks)
            else:
                raise RuntimeError("Warning: something unexpected has happened!")
        return nicks  # returns remaining nicks

    def cycle_loop(self, loop):
        """
        Produces a new sequence with a protected to un-protected 5',3' pair to start
        :param loop: List of Nucleotide of a strand that is a loop
        :return: Strand, ordered tuple of Nucleotides (Why tuple and not list?)
        """
        n0 = loop[0]
        n5 = n0
        n3 = n5.toThree
        while not (n5 in self.protected and n3 not in self.protected):
            n5 = n3
            n3 = n5.toThree
        return tuple(strandnav.getstrand(n5))

    def split2substrands(self, loop):
        """
        Consolidates the protected and unprotected regions on a strand
        Traverse a strand and returns it as an ordered list of protected and unprotected regions
        :@invariant:
          Even indices will always be protected substrands
          Odd indices will always be unprotected substrands
        :param loop: Strand, list of nucleotides where first and last elements are connected to form a cycle
        :return:
        """
        this_nucl = init_nucl = loop[0]  # find somewhere to start
        # move until we arrive at the start of a protected region in the 5'-3' direction
        while True:
            # Valid starting region is 5' no protected and 3' protected
            if this_nucl not in self.protected and this_nucl.toThree in self.protected:
                this_nucl = this_nucl.toThree
                log.dev(__name__, "Found a valid start point at {}".format(this_nucl.numid))
                break
            # Otherwise continue searching into the 3' direction
            else:
                log.dev(__name__, "{}-{} is not valid startpoint, moving in 3'".format(
                    this_nucl.numid, this_nucl.toThree.numid))
                this_nucl = this_nucl.toThree
            # Getting back to the initial nucleotide means the strand is fully protected
            #   Thus, entire strand can be immediately returned
            if this_nucl == init_nucl:
                log.dev(__name__, "split2substrands: This substrand is a fully protected loop.")
                return [loop]

        log.dev(__name__, "Starting consolidation at {}".format(this_nucl.numid))
        # From this point on, consolidate into separate strands
        init_nucl = this_nucl
        strand = []  # List of nucleotides for each strand
        cstrands = []  # List of all consolidated substrands
        while True:
            # Case 1
            # traversing a protected substrand and nucl is protected
            # add nucl to strand
            if this_nucl in self.protected and not len(cstrands) % 2:
                strand.append(this_nucl)
            # Case 2
            # traversing a unprotected substrand and nucl is unprotected
            elif this_nucl not in self.protected and len(cstrands) % 2:
                strand.append(this_nucl)
            # Case 3
            # traversing an unprotected substrand and nucl is protected
            # transition n->p: add strand to cstrands
            elif this_nucl in self.protected and len(cstrands) % 2:
                cstrands.append(strand)
                strand = [this_nucl]
            # Case 4
            # traversing a protected substrand and nucl is unprotected
            # transition p->n: add strand to cstrands
            elif this_nucl not in self.protected and not len(cstrands) % 2:
                cstrands.append(strand)
                strand = [this_nucl]
            # Loop increment step
            #   Continue in 3' direction
            this_nucl = this_nucl.toThree
            # Loop termination condition
            #   Strand loop fully traversed
            if this_nucl == init_nucl:
                cstrands.append(strand)
                log.dev(__name__, "Arrived back at {}, terminating".format(init_nucl.numid))
                break
        return cstrands

    def isfullyprotectedloop(self, loop):
        """
        Checks a strand property for being fully protected
        :param loop: List of sequential nucleotides
        :return: False as soon as unprotected nucleotide is encountered, True otherwise
        """
        nucl = loop[0]
        strand = strandnav.getstrand(nucl)
        for nucl in strand:
            if nucl not in self.protected:
                return False
        return True

    def nickstrand(self, strand, fdist, ignoremarking=False):
        """
        In a list of nucleotides, add Nicks spacing according to the input feature distance
        :param strand: List of sequential nucleotides
        :param fdist: int, feature distance, how many unedited bases should be between each Nick
        :param ignoremarking: Ignore the protected list or not
        :return: List of Nick(object) placed on the strand
        """
        nicks = []
        for n in strand:
            # Don't care about protected bases
            if ignoremarking:
                if strandnav.distfromfeature(n) >= fdist:
                    nicks.append(Nick(crossover.NuclPair(n, n.toThree)))
            # Care about protected bases
            else:
                if strandnav.distfromfeature(n) >= fdist and n not in self.protected:
                    nicks.append(Nick(crossover.NuclPair(n, n.toThree)))
        # Troubleshooting
        # log.dev(__name__, "Added {} nicks:".format(len(nicks)))
        # for n in nicks:
        #     log.dev(__name__, "obj Nick: {}".format(n.nuclpair.numid))
        return nicks

    def estimate_scaffold_count(self):
        """
        Go through the pathway for JoinedStrand 's only
        Get the first scaffold sequence
        At each JoinedStrand, subtract from the scaffold
        If scaffold length becomes negative, get the next scaffold, count the same strand again
        Continue until reaching the end of the pathway
        :return: number of scaffolds
        """
        seqnames = iter(config.AVAIL_SEQUENCES)
        seqname = next(seqnames)
        seqlen = len(sequence.getseq(seqname))
        scafcount = 1
        for s in mymath.edgelist2nodelist(self.pathway_strands):
            seqlen -= s.get_length()
            if seqlen <= 0:
                seqname = next(seqnames)
                seqlen = len(sequence.getseq(seqname))
                seqlen -= s.get_length()
                scafcount += 1

        return scafcount

    def balance(self, staple):
        """
        Takes a very short strand and attempts to bring it up to the lower limit for DNA origami (15)
        :param staple: Strand, List of nucleotides
        :return: None, performs new nicks
        """
        # Find the corresponding nick
        strand5 = strandnav.get5strand(staple)
        strand3 = strandnav.get3strand(staple)
        nick5 = self.get_nick(staple[0])
        nick3 = self.get_nick(staple[-1])

        # Undo the 5' nick
        # See if there is a better location
        for nucl in strand5:
            if strandnav.distfromfeature(nucl) >= 7:
                log.system("Valid alternative found on 5 strand")
                nick5.undo()
                self.nicks.remove(nick5)
                self.nicks.append(Nick(crossover.NuclPair(nucl, nucl.toThree)))
                return

        # If not, undo the 3' nick
        # See if there is a better location
        for nucl in strand3:
            if strandnav.distfromfeature(nucl) >= 7:
                log.system("Valid alternative found on 3 strand")
                nick3.undo()
                self.nicks.remove(nick3)
                self.nicks.append(Nick(crossover.NuclPair(nucl, nucl.toThree)))
                return

        # Raise an error
        raise RuntimeError("There were no rebalance alternatives to resolve this short strand.")

    def repair(self):
        """
        Repairs the nick map
        :return:
        """
        all_modules = self.get_modules()
        for module in all_modules:
            for nucl in module.helix.stap:
                if nucl.toThree == -1 and nucl.__strand3__.toFive == -1:
                    found = False
                    for nick in self.nicks:
                        if nucl in nick.nuclpair.to_tuple:
                            found = True
                    if not found:
                        self.nicks.append(Nick(crossover.NuclPair(nucl, nucl.__strand3__)))
        # TODO: Getting nicks by object doesn't work yet because some Nicks are not added correctly to nick_map

    def clean_merge(self):
        """
        Forces merges of short strands with adjacent strands
        :return:
        """
        all_staples = self.get_all_staples()
        for strand in all_staples:
            if len(strand) < 20:
                end5 = strandnav.goto5(strand[0])
                end3 = strandnav.goto3(strand[0])
                adj5 = end5.__strand5__
                adj3 = end3.__strand3__
                strand5 = strandnav.getstrand(adj5)
                strand3 = strandnav.getstrand(adj3)

                len5 = len3 = self.get_size()  # initiate with very high value
                if not strandnav.strandequals(strand, strand5) and not end5.is_on_gap():
                    len5 = len(strand) + len(strand5)
                elif not strandnav.strandequals(strand, strand3) and not end3.is_on_gap():
                    len3 = len(strand) + len(strand3)
                else:
                    raise RuntimeError("This short strand cannot be forcefully joined to any adjacent strand.")

                if len5 < len3:
                    dnaconnector.nuclconnect(adj5, end5)
                    log.system(
                        "[WARNING] Forcefully reconnected a short strand of length {} to an adjacent strand of length "
                        "{} for total length {}.".format(
                            len(strand), len(strand5), len(strand) + len(strand5)))
                elif len3 < len5:
                    dnaconnector.nuclconnect(end3, adj3)
                    log.system(
                        "[WARNING] Forcefully reconnected a short strand of length {} to an adjacent strand of length "
                        "{} for total length {}.".format(
                            len(strand), len(strand3), len(strand) + len(strand3)))
                else:
                    raise RuntimeError("Neither length was modified. "
                                       "This short strand cannot be forcefully joined to ady adjacent strand.")

    def clean_break(self):
        """
        Forces nicks into long strands after all other processes have finished
        :return:
        """
        pass


"""
================
Static functions
================
"""


def expand(strand):
    """
    Lengthens the strand at both endpoints to attempt to encompass more helices
    :param strand: list of Nucleotide
    :return: None, rearranges Nicks
    """
    """
    Just do this by hand for now
    """
    pass


def break_loop(refnucl):
    """
    Finds the best location to add a nick on a strand that is looped but within the bounded length of a strand
    :param refnucl:
    :return: None
    """
    if not strandnav.checkloop(refnucl):
        raise RuntimeError("Not a loop.")

    strand = strandnav.getstrand(refnucl)
    bestpos = (None, 0)  # Nucleotide, farthest distance from feature
    for nucl in strand:
        thisdist = strandnav.distfromxover(nucl)
        if thisdist > bestpos[1]:
            bestpos = (nucl, thisdist)

    if bestpos[1] <= config.LOOP_CHECK_RANGE:
        raise RuntimeError("Not finding any space on this loop, at all.")
    else:
        return Nick(crossover.NuclPair(bestpos[0], bestpos[0].toThree))


def _merge_strands(n5, n3, nicks):
    """
    Helper method to connect 2 strands and remove the nick
    :param n5: 5' Nucleotide
    :param n3: 3' Nucleotice
    :param nicks: List of Nick(object)
    :return: Updated list
    """
    for nick in nicks:  # remove the nick from list
        if nick.nuclpair.to_tuple == (n5, n3):
            nick.undo()
            nicks.remove(nick)
            break
    return nicks


def nickstostrands(nicks):
    """
    Converts a list of nicks to a list of strands
    :param nicks: List of Nick(object)
    :return: List of strands (list of Nucleotide) derived from list of Nick
    """

    strands = mymath.FakeSet([])
    if not nicks:
        # log.debug(__name__, nicks)
        return strands
    for nick in nicks:
        strands.add(tuple(strandnav.getstrand(nick.nuclpair.n5)))
        strands.add(tuple(strandnav.getstrand(nick.nuclpair.n3)))
    return strands


def excise(strands):
    """
    Splits a list of substrands of a long strand into multiple distinct strands
    Produces a list of nick objects from a list of ordered substrands
    @invariant:
    all strands are already at a protected/unprotected transition
    this operation is only safe to do if that is True
    :param strands: List of substrands consolidated from same original strand.
        Each substrand is a list of nucleotides.
    :return: List of Nick
    """
    nicks = []
    for s in strands:
        init = s[0]
        n5 = init.toFive
        n3 = init
        nicks.append(Nick(crossover.NuclPair(n5, n3)))
    return nicks
