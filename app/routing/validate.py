"""
Created by Daniel Fu (Reif Lab, Duke University) at 6/27/2020

Name        : validate.py
Project     : CADAxiSDNA
Description : Methods that run at the end of the origami generation to make sure everything was OK
Interpreter : Python 3.7.4
"""

from app import config
from .helper import strandnav, mymath, log
from . import nucleotide, sequence, modules, path
from .exceptions import *
import numpy as np
import itertools


def setmsg(name, description, value):
    """
    Small helper function to save how each validation is displayed
    :param name: Name of what was validated
    :param description: A description of the process
    :param value: Truth value
    :return: Format the output
    """
    if value:
        value = "PASS"
    else:
        value = "FAIL"
    return {"Name": name, "Description": description, "Value": value}

class Mixin:

    def validate(self):
        log.system("Running integrity checks.")

        validations = {}

        validations['sequences'] = setmsg("Applied sequences match input",
                                          "Validate that scaffold sequences were correctly applied.",
                                          self.validate_sequences())

        if config.ROUTING.lower().startswith('modular'):  # Only if origami is partitioned
            validations['partitions'] = setmsg("Partition count matches scaffold nicks",
                                               "Validate that the number of modular origami segments (determined by "
                                               "which scaffold crossovers were skipped) matches the number of scaffold "
                                               "start positions (determined by the number of nick points in the "
                                               "scaffold strands.",
                                               self.validate_partitions())

        validations['connections'] = setmsg("No strand continuity problems",
                                            "Validate that all nucleotides on each strand are all continuously "
                                            "connected in both directions 5'-3' and 3'-5'. Queries strands from "
                                            "helices",
                                            self.validate_strand_continuity)

        validations['balanced_nicking'] = setmsg("Nicks are balanced",
                                                 "Validate that the strand traversing toThree or toFive attributes of "
                                                 "each nucleotide on both sides of a nick is -1 (indicating no "
                                                 "connection.",
                                                 self.validate_balanced_nicking())
        # Not implemented
        validations['antiparallel'] = setmsg("Adjacent helices anti-parallel",
                                             "Validate that adjacent helices are anti-parallel",
                                             self.validate_antiparallel())

        validations['strand_query'] = setmsg("Strand query consistency",
                                             "Validate that strands queried from the helix property of modules "
                                             "matches the strands queried from looking for 5' nicks.",
                                             self.validate_strand_queries())

        validations['dupl_nucl'] = setmsg("No duplicate nucleotides",
                                          "Validate that no strands share the same nucleotides",
                                          self.validate_dupl_strand_nucl() and self.validate_dupl_helix_nucl())

        validations['size'] = setmsg("Size check",
                                     "Validate that the size of the structure matches when queried from strands, "
                                     "nucleotides, and modules.",
                                     self.validate_bp_total())

        validations['unassigned_nucl'] = setmsg("All nucleobases labelled",
                                                "Checks that all nucleotides have a valid base letter.",
                                                self.validate_unassigned_nucleotides())

        validations['nucl_complement'] = setmsg("All strands complementary",
                                                "Validates that all nucleotides have a set complementary nucleotide.",
                                                self.validate_complements())

        validations['xover_distance'] = setmsg("No stretched crossovers",
                                               "Validates that all crossovers have reasonable bond lengths measured "
                                               "in Euclidean distance.",
                                               self.validate_xover_distance())

        validations['nucl_topology'] = setmsg("Nucleotide-strand topology",
                                              "Every nucleotide on the strand is topologically linked to that strand",
                                              self.validate_nucl_top())

        validations['short_regions'] = setmsg("Sufficient region lengths",
                                              "Regions (the uninterrupted length of a double helix with no nicks or "
                                              "crossovers) that are too short are avoided.",
                                              self.validate_short_regions())

        validations['short_strands'] = setmsg("Sufficient staple lengths",
                                              "Staple lengths should not be lower than the config lower bound",
                                              self.validate_short_staples())

        validations['long_oligos'] = setmsg("Within IDT upper bound (90)",
                                             "Can be synthesized by IDT.",
                                             self.validate_staple_high_bound())

        validations['short_oligos'] = setmsg("Within IDT lower bound (15)",
                                             "Can be synthesized by IDT.",
                                             self.validate_staple_low_bound())

        validations['scloops'] = setmsg("Loops in scaffold routing",
                                        "Validates that scaffold crossovers are placed correctly",
                                        self.validate_scloops())

        validations['stloops'] = setmsg("Loops in staple routing",
                                        "Validates that staple nicking proceeded correctly",
                                        self.validate_stloops())

        validations['loops'] = setmsg("Loops in either routing",
                                      "Validates that all routings proceeded correctly",
                                      self.validate_loops())

        validations['scaffold_anti_parallel'] = setmsg("All scaffold crossovers anti-parallel.",
                                                       "Validates that the scaffold routing only uses anti-parallel "
                                                       "crossovers. Mixing with parallel crossovers will create "
                                                       "isolated loops.",
                                                       self.validate_scaf_antiparallel())

        validations['scaffold_xovers'] = setmsg("Scaffold crossovers match removed list.",
                                                "Checks for crossovers between modules and compares it to crossovers "
                                                "that were supposed to be skipped.",
                                                self.validate_scaf_xovers())

        validations['xover_overlaps'] = setmsg("No crossover overlaps.",
                                               "Validates that no crossovers were formed on the same nucleotides such "
                                               "that it effectively skips a strand.",
                                               self.validate_xover_overlaps())

        validations['single_scaf_xover_path'] = setmsg("Scaffold crossovers per pair are unique.",
                                                       "Validates that there is not more than 1 crossover between any "
                                                       "pair of nearest neighbor helices.",
                                                       self.validate_single_scaf_xover_path())

        validations['gap_staple_terminates'] = setmsg("Gap staples are nicked.",
                                                       "Gap staples should be nicked in acyclic structures otherwise "
                                                       "they are likely to form poor length regions.",
                                                      self.validate_gap_staple_terminates())

        validations['nick_map'] = setmsg("Are all nicked locations saved.",
                                         "Check whether nick_map matches all placed nick.",
                                         self.validate_nick_map())

        validations['xovers_complete'] = setmsg("Are all crossovers saved.",
                                                "Checks whether objxovers matches all placed crossovers.",
                                                self.validate_xovers_complete())

        validations['seams'] = setmsg("Seams detected.",
                                      "Checks for scaffold and staple crossovers too close together on a common helix.",
                                      self.validate_seams())

        validations['sparse_xovers'] = setmsg("No sparsely connected helices",
                                              "Checks whether all nearest-neighbor connected helices have at least one "
                                              "staple crossover",
                                              self.validate_sparse_xovers())

        flag = True

        for key in validations:
            log.system("{:.<50}{}".format(validations[key]['Name'], validations[key]['Value']))
            if validations[key]['Value'] != 'PASS':
                flag = False
        if flag:
            log.system("All OK.")
        else:
            log.log_warning("Problems detected. Check log.")

    def validate_sparse_xovers(self):
        """
        Check if any pair of adjacent modules have too few crossovers (only a scaffold)
        :return: True/False
        """
        sparse_flag = False
        for key in self.man_connections:
            module1 = self.planes[key[0]].modules[key[1]]
            for val in self.man_connections[key]:
                module2 = self.planes[val[0]].modules[val[1]]
                js1 = module1.get_top_strand()
                js2 = module2.get_top_strand()
                # If there is not any FullXover of strandtype staple spanning the two JoinedStrands, the it is spare
                if not any(xobj.spans(js1, js2) and xobj.strandtype == 'stap' for xobj in module1.objxovers):
                    sparse_flag = True
                    log.system("[WARNING]: Crossovers between {} and {} are sparse".format(module1, module2))
        if sparse_flag:
            return False
        else:
            return True

    def validate_seams(self):
        """
        Checks if there are seams throughout the structure
        A seam is a staple and scaffold crossover that share a common helix and are spaced too close together
        :return: True/False
            Also prints out cases that are too close together
        """
        seam_flag = True
        all_objxovers = self.get_all_xovers()
        for xover in all_objxovers:
            if not xover.strandtype == "scaf":
                continue
            (n1, n2) = (xover.n1, xover.n2)
            (n1_5, n1_3) = n1.to_tuple
            (n2_5, n2_3) = n2.to_tuple

            for i in range(3):
                n1_5c = n1_5.Comp
                n1_3c = n1_3.Comp
                n2_5c = n2_5.Comp
                n2_3c = n2_3.Comp
                check = (strandnav.isxover(n1_5c), strandnav.isxover(n1_3c),
                         strandnav.isxover(n2_5c), strandnav.isxover(n2_3c))
                if any(check):
                    log.system("A seam is at the following crossover between modules {} and {}: {}".format(n1.up(),
                                                                                                           n2.up(),
                                                                                                           xover))
                    seam_flag = False
                n1_5 = n1_5.__strand5__
                n1_3 = n1_3.__strand3__
                n2_5 = n2_5.__strand5__
                n2_3 = n2_3.__strand3__

        if seam_flag:
            return True
        else:
            return False

    def validate_xovers_complete(self):
        """
        Checks if all crossovers detected by parsing through the structure are the same that have been logged by
        all modules
        :return: True/False
        """
        cnt_nuclnotfound = 0
        cnt_danglingxover = 0
        all_objxovers = self.get_all_xovers()
        all_strands = self.get_all_strands()
        for strand in all_strands:
            for nucl in strand:
                if strandnav.isxover(nucl):
                    try:
                        self.findxover(nucl)
                    except NuclNotFound:
                        cnt_nuclnotfound += 1
                        log.system("[WARNING]: A detected crossover was not listed on any module.")
        for xover in all_objxovers:
            (n1, n2) = (xover.n1, xover.n2)
            (n1_5, n1_3) = n1.to_tuple
            (n2_5, n2_3) = n2.to_tuple
            if not (n1_5.toThree == n2_3 and n1_3.toFive == n2_5):
                log.system("[WARNING]: A listed crossover is not applied.")
                cnt_danglingxover += 1

        if cnt_danglingxover>0 or cnt_nuclnotfound>0:
            return False
        else:
            return True

    def validate_nick_map(self):
        # print("All the nicks:", self.nicks)

        all_modules = self.get_modules()
        all_toThree_nicks = []
        for m in all_modules:
            for nucl in m.helix.stap:
                if nucl.toThree == -1:
                    all_toThree_nicks.append(nucl)

        if len(self.nicks) == len(all_toThree_nicks):
            return True
        else:
            log.system("There are {} saved nicks but {} detected nicks.".format(len(self.nicks), len(all_toThree_nicks)))

    def validate_short_regions(self):
        """
        Checks for short regions
        :return: True/False
        """
        all_regions = self.get_all_regions()
        for region in all_regions:
            if len(region) < 3:
                return False
        return True

    def validate_short_staples(self):
        """
        Checks for short staples
        :return: True/False
        """
        all_staples = self.get_all_staples()
        for staple in all_staples:
            if len(staple) < config.LENLOW:
                return False
        return True

    def validate_gap_staple_terminates(self):
        """
        Checks whether staple strands are terminating (nick) at gap
        :return: True/False
        """
        all_modules = self.get_modules()
        for module in all_modules:
            if type(module) is modules.ArcModule:
                endpt1 = module.helix.stap[0]
                endpt2 = module.helix.stap[-1]
                if (endpt1.toThree != -1 and endpt1.toFive != -1) or (endpt2.toThree != -1 and endpt2.toFive != -1):
                    return False
        return True

    def validate_scloops(self):
        """
        Checks if there are still loops in the scaffold routing
        Most likely implies problems with scaffold crossover placement
        :return: True/False
        """
        scaf_loops, _ = self.get_all_strand_loops(strandtype='scaf')
        # log.system("Content of only scaf_loops", scaf_loops)
        if scaf_loops:
            return False
        else:
            return True

    def validate_stloops(self):
        """
        Checks if there are still loops in the staple routing
        Most likely implies problems with the nicking algorithm
        :return: True/False
        """
        stap_loops, _ = self.get_all_strand_loops(strandtype='stap')
        # log.system("Content of only stap_loops", stap_loops)
        if stap_loops:
            return False
        else:
            return True

    def validate_loops(self):
        """
        Checks if there are still loops in the routing
        :return: True/False
        """
        stap_loops, _ = self.get_all_strand_loops(strandtype='stap')
        scaf_loops, _ = self.get_all_strand_loops(strandtype='scaf')
        # log.system("Content of both scaf_loops", scaf_loops)
        # log.system("Content of both stap_loops", stap_loops)
        if stap_loops or scaf_loops:
            return False
        else:
            return True

    def validate_staple_high_bound(self):
        """
        Checks if all staples are under the custom oligo length upper limit (90)
        :return: True/False
        """
        all_staples = self.get_all_staples()
        for strand in all_staples:
            if len(strand) > 90:
                return False
        return True

    def validate_staple_low_bound(self):
        """
        Checks if all staples are above the custom oligo length lower limit (15)
        :return: True/False
        """
        all_staples = self.get_all_staples()
        for strand in all_staples:
            if len(strand) < 15:
                # 5', 3' ends of strand
                sn5 = strandnav.goto5(strand)
                sn3 = strandnav.goto3(strand)
                try:
                    adj5 = strandnav.getstrand(sn5.__strand5__)
                except AttributeError:
                    adj5 = []
                try:
                    adj3 = strandnav.getstrand(sn3.__strand3__)
                except AttributeError:
                    adj3 = []
                log.system("[ERROR]: Critically short strands detected: "
                           "\n\t5: {}"
                           "\n\ts: {}"
                           "\n\t3: {}".format(len(adj5),
                                                 len(strand),
                                                 len(adj3)))
                return False
        return True

    def validate_scaf_xovers(self):
        """
        Checks if each adjacent scaffold strand pair has a crossover
        (Is the scaffold strand fully routed through the structure?)
        :return: True/False
        """
        log.debug("Validating scaffold routing via crossovers.")
        for edge in self.pathway:
            mdx1 = edge[0]
            md1 = self.planes[mdx1[0]].modules[mdx1[1]]
            if edge[1]:
                mdx2 = edge[1]
                md2 = self.planes[mdx2[0]].modules[mdx2[1]]
            else:
                continue
            found_flag = False
            for nucl in md1.helix.scaf:
                try:
                    if nucl.toFive.up().up() == md2:
                        found_flag = True
                        break
                except AttributeError:  # Ran into a break first
                    pass
            if found_flag:
                log.debug("{:40} {} and {}".format("OK yes there's a crossover between", md1, md2))
            else:
                if path.Edge(md1, md2) not in self.removed.getsorted():
                    log.debug("A skipped scaffold crossover does not match the removed list: {}".format(
                        self.removed.getsorted()))
                    return False
                log.debug("{:40} {} and {}".format("No crossover between", md1, md2))
        return True

    def validate_xover_overlaps(self):
        """
        Checks if any crossovers overlap (more than one crossover uses the same nucleotide)
        :return: True/False
        """
        all_modules = self.get_modules()
        all_xovers = mymath.FakeSet([])
        for module in all_modules:
            [all_xovers.add(xobj) for xobj in module.objxovers]
        all_nucls = []
        for x in all_xovers:
            [all_nucls.append(n.numid) for n in [x.n1.n5, x.n1.n3, x.n2.n5, x.n2.n3]]
        if len(all_nucls) != len(set(all_nucls)):
            return False
        else:
            return True

    def validate_single_scaf_xover_path(self):
        """
        Checks if there are more than one scaffold crossover between helices
        :return: True/False
        """
        for jstrand in self.strands:
            xover_counter = 0
            init_nucl = jstrand.modules[0].helix.scaf[0]
            nucl = init_nucl
            while True:
                try:
                    if nucl.toThree.get_top_strand() != jstrand:
                        xover_counter += 1
                except AttributeError:
                    pass
                nucl = nucl.__strand3__
                if nucl == init_nucl:
                    break
            if xover_counter > 2:
                return False
        return True

    def validate_scaf_antiparallel(self):
        """
        Checks that all scaffold crossovers are anti-parallel
        :return: True/False
        """
        all_modules = self.get_modules()
        all_xovers = []
        for m in all_modules:  # Get all relevant crossovers
            for x in m.objxovers:
                if x.strandtype == 'scaf':
                    all_xovers.append(x)
        while True:  # Clear duplicates
            for x1, x2 in itertools.permutations(all_xovers, 2):
                if x1 == x2:
                    all_xovers.remove(x2)
                    break
            if all(not x1 == x2 for x1, x2 in itertools.permutations(all_xovers, 2)):
                break
        if all(x.is_anti_parallel() for x in all_xovers):  # Check if anti-parallel
            return True
        else:
            return False

    def validate_balanced_nicking(self):
        """
        Check that every nick has a toThree = -1 and toFive = -1 on both sides
        :return: True/False
        """
        all_nucl = self.get_all_helix_nucl()
        for nucl in all_nucl:
            if nucl.toThree == -1:
                if nucl.__strand3__.toFive != -1:
                    return False
            if nucl.toFive == -1:
                if nucl.__strand5__.toThree != -1:
                    return False
        return True

    def validate_nucl_top(self):
        """
        Check that every nucl on the strand calls the strand
        :return: True/False
        """
        staple_strands = self.get_all_staples()
        for s in staple_strands:
            for nucl in s:
                cstrand = nucl.up().stap
                if nucl not in cstrand:
                    return False
        return True

        # Validate module connections
        # Shows the connection angle of adjacent modules
        # debug.validate_module_connections(dnaorigami, 30)

        # Validate basis direction/consistency
        # Check for adjacent nucleotides that have different right-handed twist directions
        # for strand in dnaorigami.get_all_strands():
        #     refnucl = strand[0]
        #     debug.validate_basis_direction(refnucl)

        # Validate twist for selected modules
        # (a, b) = (0, 0)
        # refnucl = dnaorigami.planes[a].modules[b].helix.stap[0]
        # debug.validate_twist(refnucl)
        #
        # (a, b) = (0, 8)
        # refnucl = dnaorigami.planes[a].modules[b].helix.stap[0]
        # debug.validate_twist(refnucl)
        #
        # (a, b) = (1, 0)
        # refnucl = dnaorigami.planes[a].modules[b].helix.stap[0]
        # debug.validate_twist(refnucl)
        #
        # (a, b) = (1, 8)
        # refnucl = dnaorigami.planes[a].modules[b].helix.stap[0]
        # debug.validate_twist(refnucl)
        # Check if there are crossovers with NuclPair theta that are actually not aligned due to mismatched basis vectors
        # dnaorigami.validate_xover_basis()

    def validate_dupl_strand_nucl(self):
        """
        Integrity check for same nucleotide listed twice in strand topology
        :return: print duplicate positions as numid Nucleotide and False, True if None
        """
        all_nucl = self.get_all_strand_nucl()
        for n in all_nucl:
            if all_nucl.count(n) > 1:
                log.system("{} is observed {} times in strands.".format(n.numid, all_nucl.count(n)))
                return False
        return True

    def validate_dupl_helix_nucl(self):
        """
        Integrity check for same nucleotide listed twice in strand topology
        :return: print duplicate positions as numid Nucleotide and False, True if None
        """
        all_nucl = self.get_all_helix_nucl()
        for n in all_nucl:
            if all_nucl.count(n) > 1:
                log.system("{} is observed {} times in strands.".format(n.numid, all_nucl.count(n)))
                return False
        return True

    def validate_partitions(self):
        """
        Compares that the number of scaffold nicks (currently added manually) matches the number of partitions formed
        from crossover edges that were removed (as was determined automatically due to automatic placement of scaffold
        sequences)
        :return: True if match, False if not matching
        """
        # Get removed list from break plane procedure.")
        starts = self.get_scaf_starts()
        removed = self.break_plane_asym()
        partitions = self.partition_from_removed(removed)

        if len(partitions) == len(starts):
            return True
        return False


    def validate_antiparallel(self):
        """
        Checks that all anti-parallel helices have opposing 5'-3' directions
        :return:
        """
        return True

    def validate_strand_queries(self):
        """
        Integrity check for strand topology methods.
        Checks if strands parsed with query.get_all_strands() contain the same nucleotides as
            those specified in scaf and stap lists of DNAHelix object
        :return: Prints a list of differing nucleotides. Empty is nothing is wrong.
        """
        helix_nucl = self.get_all_helix_nucl()
        strands = self.get_all_strands()
        strand_nucl = []
        for s in strands:
            for nucl in s:
                strand_nucl.append(nucl)
        log.system("Helix Nucls: {} , Strand Nucls: {}".format(len(helix_nucl), len(strand_nucl)))
        # in helix, not strand
        h_not_s = []
        for nucl in helix_nucl:
            if nucl not in strand_nucl:
                h_not_s.append(nucl)
        # in strand, not helix
        s_not_h = []
        for nucl in strand_nucl:
            if nucl not in helix_nucl:
                s_not_h.append(nucl)

        if h_not_s or s_not_h:
            for h in h_not_s:
                log.system("Found {} in helix but not in strand".format(h.numid))
            for s in s_not_h:
                log.system("Found {} in strand but not in helix".format(s.numid))
            flag = False

        else:  # not h_not_s and not s_not_h:
            log.system("Helix and strand nucleotides are consistent.")
            flag = True

        return flag

    def validate_sequences(self):
        """
        Checks that the scaffold strands matches the set sequences
        :return: True if match, False if error
        """
        # Get the scaffold starts
        starts = self.get_scaf_starts()
        num_starts = len(starts)

        # Iterate through each scaffold
        for seq_num in range(num_starts):
            seqname = config.AVAIL_SEQUENCES[seq_num]
            this_seq = sequence.getseq(seqname)
            this_strand = strandnav.getstrand(starts[seq_num])
            len_strand = len(this_strand)
            for nuclidx in range(len_strand):
                if this_seq[nuclidx] != this_strand[nuclidx].nucl:
                    return False
        return True

    def validate_strand_continuity(self):
        """
        Integrity check after edits such as crossovers, insertions, and deletions
        Checks every nucleotide to confirm that its 5', 3', and complementary neighbors all link back to itself
        This ensures there are no problems with strand or ring topology
        :return: True if nothing is wrong. Raises error otherwise
        """
        all_nucl = self.get_all_helix_nucl()

        # check if every connection is bidirectional
        for nucl in all_nucl:
            try:
                if nucl != nucl.toThree.toFive:
                    log.system("[DISCONTINUOUS] Numid: {}, Printid: {}".format(nucl.numid, nucl.printid))
            except:
                log.system("Numid: {}, Print: {} toThree is a nick.".format(nucl.numid, nucl.printid))

            try:
                if nucl != nucl.toFive.toThree:
                    log.system("[DISCONTINUOUS] Numid: {}, Printid: {}".format(nucl.numid, nucl.printid))
            except:
                log.system("Numid: {}, Print: {} toFive is a nick.".format(nucl.numid, nucl.printid))
            try:
                if nucl != nucl.Comp.Comp:
                    log.system("[BAD HYBRIDIZATION] Numid: {}, Printid: {}".format(nucl.numid, nucl.printid))
            except:
                log.system("Numid: {}, Print: {} Comp missing".format(nucl.numid, nucl.printid))
        log.system("All good!")
        return True

    def validate_bp_total(self):
        """
        Compare the baseIDtotal according to strand lengths, nucleotide count, and module.bp
        """
        strand_baseIDtotal = 0
        for s in self.get_all_strands():
            strand_baseIDtotal = strand_baseIDtotal + len(s)
        nucl_baseID_total = len(self.get_all_helix_nucl())
        module_baseID_total = sum([m.bp for m in self.get_modules()]) * 2
        # log.system("[Totals] Strand: {}, Helix Nucl: {}, Module Bp: {}".format(
        #     strand_baseIDtotal, nucl_baseID_total, module_baseID_total))
        return strand_baseIDtotal == nucl_baseID_total == module_baseID_total

    # def validate_insertions(self):
    #     """
    #     UNUSED
    #     Validate number of insertions
    #     :return:
    #     """
    #     log.system("\t# Validate number of insertions")
    #     all_nucl = self.get_all_helix_nucl()
    #     num_insertions = 0
    #     for nucl in all_nucl:
    #         if nucl.numid > 5000:
    #             num_insertions += 1
    #             log.system("Insertion at {} / {}".format(nucl.numid, nucl.Comp.numid))
    #     log.system("Number of insertions = {}".format(num_insertions))

    # ===========
    # sequence.py
    # ===========
    def validate_unassigned_nucleotides(self):
        """
        Checks for unassigned nucleotides
        :return: Return True if none unassigned, False otherwise and print numid
        """
        for nucl in self.get_all_helix_nucl():
            if nucl.nucl == 'X':
                log.system(nucl.numid)
                return False
        return True

    def validate_complements(self):
        """
        Check if the sequence and its complement were applied correctly
        :return:
        """
        for strand in self.get_all_strands():
            for nucl in strand:
                comp = nucl.Comp
                nnucl = nucl.nucl
                ncomp = comp.nucl
                flag = nnucl == nucleotide.comp_nucl(ncomp)
                # log.system("{} {} {}".format(nnucl, ncomp, flag))
                if not flag:
                    return False
        return True

    """
    Strand Debugging
    """

    # def validate_module_connections(self, reference=None):
    #     """
    #     Gets the endpoints of all the modules and checks
    #     Run this for structures where a single strand is composed of multiple modules
    #     Run before applying crossovers and after applying module connections
    #     :param dnaorigami: Origami(object)
    #     :return: None, prints results in console output only
    #     """
    #     # Check if each dtheta equals the desired reference
    #     if reference:
    #         value = reference
    #     for plane in self.planes:
    #         for module in plane.modules:
    #             # Get the endpoints of the module
    #             end1 = module.helix.scaf[0].node
    #             end2 = module.helix.scaf[-1].node
    #             scaf5 = end1.scaf
    #             scaf3 = end2.scaf
    #             stap5 = end2.stap
    #             stap3 = end1.stap
    #
    #             # If the strand is connected (and not terminal)
    #             if scaf5.toFive != -1:
    #                 dt = mytrig.angle_diff(scaf5.toFive.theta, scaf5.theta)
    #                 nucl = scaf5
    #                 nucl5 = nucl.toFive
    #                 nucl55 = nucl.toFive.toFive
    #                 nucl3 = nucl.toThree
    #                 nucl33 = nucl.toThree.toThree
    #                 log.system("scaf5: numid: {}, dTheta: {}\t\t\t{}".format(scaf5.numid, dt, dt == value))
    #                 log.system("\t{:5.1f} | {:5.1f} | {:5.1f} | {:5.1f} | {:5.1f}".format(
    #                     nucl55.theta, nucl5.theta, nucl.theta, nucl3.theta, nucl33.theta))
    #             if scaf3.toThree != -1:
    #                 dt = mytrig.angle_diff(scaf3.theta, scaf3.toThree.theta)
    #                 nucl = scaf3
    #                 nucl5 = nucl.toFive
    #                 nucl55 = nucl.toFive.toFive
    #                 nucl3 = nucl.toThree
    #                 nucl33 = nucl.toThree.toThree
    #                 log.system("scaf3: numid: {}, dTheta: {}\t\t\t{}".format(scaf3.numid, dt, dt == value))
    #                 log.system("\t{:5.1f} | {:5.1f} | {:5.1f} | {:5.1f} | {:5.1f}".format(
    #                     nucl55.theta, nucl5.theta, nucl.theta, nucl3.theta, nucl33.theta))
    #             if stap5.toFive != -1:
    #                 dt = mytrig.angle_diff(stap5.toFive.theta, stap5.theta)
    #                 nucl = stap5
    #                 nucl5 = nucl.toFive
    #                 nucl55 = nucl.toFive.toFive
    #                 nucl3 = nucl.toThree
    #                 nucl33 = nucl.toThree.toThree
    #                 log.system("stap5: numid: {}, dTheta: {}\t\t\t{}".format(stap5.numid, dt, dt == value))
    #                 log.system("\t{:5.1f} | {:5.1f} | {:5.1f} | {:5.1f} | {:5.1f}".format(
    #                     nucl55.theta, nucl5.theta, nucl.theta, nucl3.theta, nucl33.theta))
    #             if stap3.toThree != -1:
    #                 dt = mytrig.angle_diff(stap3.theta, stap3.toThree.theta)
    #                 nucl = stap3
    #                 nucl5 = nucl.toFive
    #                 nucl55 = nucl.toFive.toFive
    #                 nucl3 = nucl.toThree
    #                 nucl33 = nucl.toThree.toThree
    #                 log.system("stap3: numid: {}, dTheta: {}\t\t\t{}".format(stap3.numid, dt, dt == value))
    #                 log.system("\t{:5.1f} | {:5.1f} | {:5.1f} | {:5.1f} | {:5.1f}".format(
    #                     nucl55.theta, nucl5.theta, nucl.theta, nucl3.theta, nucl33.theta))
    #
    # def validate_twist(refnucl):
    #     """
    #     Gets the strand from a reference nucleotide
    #     Checks for adjacent nucleotides that are offset by too large of an angle
    #     :param refnucl: object Nucleotide
    #     :return: True if no errors, False otherwise
    #     """
    #     count = 0
    #     strand = strandnav.absgetstrand(refnucl)
    #     for nucl in strand:
    #         nucl3 = nucl.__strand3__
    #         try:
    #             t3 = mytrig.acute_diff(nucl3.theta, nucl.theta)
    #         except AttributeError:
    #             pass
    #         if t3 > 40:
    #             if nucl.up() == nucl3.up():
    #                 flag = True
    #             else:
    #                 flag = False
    #             log.system("Warning  ({:6}) @ {:5} ({:6}) {}|{:5} ({:6}) {} > {}".format(
    #                 t3, nucl.numid, nucl.theta, nucl.up(), nucl3.numid, nucl3.theta, nucl3.up(), flag))
    #             count += 1
    #     return count
    #
    # def validate_basis_direction(refnucl):
    #     """
    #     Gets the strand from a reference nucleotide
    #     Checks for inconsistent basis vectors and inconsistent theta direction
    #     :param refnucl: reference Nucleotide(object)
    #     :return: True if no errors, False if errors and prints to console
    #     """
    #     strand = strandnav.absgetstrand(refnucl)
    #     for nucl in strand:
    #         try:
    #             # Get a range of 5 nucleotides
    #             nucl3 = nucl.__strand3__
    #             nucl33 = nucl3.__strand3__
    #             nucl5 = nucl.__strand5__
    #             nucl55 = nucl5.__strand5__
    #             t3 = mytrig.angle_diff(nucl.theta, nucl3.theta)
    #             t33 = mytrig.angle_diff(nucl3.theta, nucl33.theta)
    #             t5 = mytrig.angle_diff(nucl5.theta, nucl.theta)
    #             t55 = mytrig.angle_diff(nucl55.theta, nucl5.theta)
    #             # Basis vectors of the middle 3
    #             b5 = np.linalg.norm(nucl5.bX - nucl.bX)
    #             b3 = np.linalg.norm(nucl3.bX - nucl.bX)
    #             if t5 * t3 < 0:  # If theta is not monoincreasing
    #                 log.system(
    #                     "          {:5} | {:^5} | {:^5} | {:^5} | {:^5} | {:^5}".format("numid", "5''", "5'", "0", "3'",
    #                                                                                     "3''"))
    #                 if mymath.basis_dist(nucl5.bX, nucl.bX, 1):
    #                     nt2 = nucl5.theta
    #                     nt1 = mymath.lh2rh2(nucl.theta)
    #                     diffnt = abs(mytrig.angle_diff(nt2, nt1))
    #                 elif mymath.basis_dist(nucl3.bX, nucl.bX, 1):
    #                     nt2 = nucl3.theta
    #                     nt1 = mymath.lh2rh2(nucl.theta)
    #                     diffnt = abs(mytrig.angle_diff(nt2, nt1))
    #                 else:
    #                     raise RuntimeError("Something else went wrong.")
    #                 log.system("Warning @ {:5} | {:5.1f} | {:5.1f} | {:5.1f} | {:5.1f} | {:5.1f} | {:5.1f}".format(
    #                     nucl.numid, nucl55.theta, nucl5.theta, nucl.theta, nucl3.theta, nucl33.theta, diffnt))
    #                 log.system("{:10} {:40} {:10} {:40} {!r:6} {!r:6}".format("5", str(nucl5.bX), str(nucl5.bY),
    #                                                                      str(nucl.duvec),
    #                                                                      b5 > 1, nucl5.direction))
    #                 log.system(
    #                     "{:10} {:40} {:10} {:40} {!r:6} {!r:6}".format("0", str(nucl.bX), str(nucl.bY), str(nucl.duvec),
    #                                                                    None, nucl.direction))
    #                 log.system("{:10} {:40} {:10} {:40} {!r:6} {!r:6}".format("3", str(nucl3.bX), str(nucl3.bY),
    #                                                                      str(nucl.duvec),
    #                                                                      b3 > 1, nucl3.direction))
    #         except AttributeError:
    #             pass

    # =============================================================================
    # crossover.py
    # =============================================================================
    # def validate_theta_vector(self):
    #     """
    #     While valid crossover positions along a helix are determined by only their theta value, that value may not be
    #     accurate due to inconsistencies with the nucleotide's basis vectors.
    #     Compare:
    #         1: The vector as given by calculating nucl_vec from the theta and basis vectors of the Nucleotide(object)
    #         2: The vector as given by the center to center distance between the two NuclPairs of the FullXover
    #     Convert each of these to spherical coordinate and compare their Phi. Check the difference between these two Phi
    #     angles. If the difference is large, it means that the chosen valid crossover is probably at an incorrect
    #     rotation angle along the strand.
    #     :return: print out results to console
    #     """
    #     for plane in self.planes:
    #         for module in plane.modules:
    #             for xobj in module.objxovers:
    #                 nuclpair1 = xobj.n1
    #                 nuclpair2 = xobj.n2
    #                 # 1: nucl_vec
    #                 np1v = mymath.unitvector(nuclpair1.vector)
    #                 np2v = mymath.unitvector(nuclpair2.vector)
    #                 np1v_sph = mymath.cart2sph(*np1v)
    #                 np2v_sph = mymath.cart2sph(*np2v)
    #                 np1vphi = np.degrees(np1v_sph[2])
    #                 np2vphi = np.degrees(np2v_sph[2])
    #                 # 2: center to center
    #                 np1c = nuclpair1.coords
    #                 np2c = nuclpair2.coords
    #                 np1cv = mymath.unitvector(np2c - np1c)
    #                 np2cv = mymath.unitvector(np1c - np2c)
    #                 np1c_sph = mymath.cart2sph(*np1cv)
    #                 np2c_sph = mymath.cart2sph(*np2cv)
    #                 np1cphi = np.degrees(np1c_sph[2])
    #                 np2cphi = np.degrees(np2c_sph[2])
    #                 # Compare
    #                 np1d = mytrig.angle_diff(np1cphi, np1vphi)
    #                 np2d = mytrig.angle_diff(np2cphi, np2vphi)
    #
    #                 # Print results
    #                 n1module = str(nuclpair1.n5.up().up())
    #                 log.system("NuclPair1 {} on {}".format(nuclpair1.numid, n1module))
    #                 log.system("NuclVec: {:4.2f}, CenterVec: {:4.2f}, Difference: {:4.2f}".format(np1vphi, np1cphi, np1d))
    #                 n2module = str(nuclpair2.n5.up().up())
    #                 log.system("NuclPair2 {} on {}".format(nuclpair2.numid, n2module))
    #                 log.system("NuclVec: {:4.2f}, CenterVec: {:4.2f}, Difference: {:4.2f}\n".format(np2vphi, np2cphi, np2d))
    #
    # def validate_nucl_theta(self):
    #     """
    #     Checks that the theta attribute of each Nucleotide(object) can be found in reverse from the vector
    #     To create a Nucleotide, input theta is converted to a nucl_vec using the basis vectors bX, bY with normal duvec
    #     If theta solved from nucl_vec is equal to input theta then return True
    #     :return: True if match and continue for all Nucleotides, False otherwise and printout error
    #     """
    #     for nucl in self.get_all_helix_nucl():
    #         in_theta = nucl.theta
    #         nucl_vec = nucl.nucl_vec
    #         (bx, bY) = (nucl.bX, nucl.bY)
    #         pass  # remains unfinished
    #
    # def validate_xover_basis(self):
    #     """
    #     Checks whether crossovers are created between strands of the same basis
    #     And if not, how far away their crossover theta are if on the same theta
    #     :return: printout
    #     """
    #     count = 0
    #     for plane in self.planes:
    #         for module in plane.modules:
    #             for xover in module.objxovers:
    #                 n1 = xover.n1
    #                 n1n5 = n1.n5
    #                 n1n3 = n1.n3
    #                 n2 = xover.n2
    #                 n2n5 = n2.n5
    #                 n2n3 = n2.n3
    #                 if mymath.basis_dist(n1n5.bX, n2n3.bX, 1) or mymath.basis_dist(n1n3.bX, n2n5.bX, 1):
    #                     log.system("Crossover confirmed to having unmatching basis vectors.")
    #                     log.system("{:15} | {:15}".format("Native basis", "Swapped basis"))
    #                     log.system("{:15} | {:15}".format(round(n1.theta, 1), round(mymath.lh2rh2(n1.theta), 1)))
    #                     log.system("{:15} | {:15}".format(round(n2.theta, 1), round(mymath.lh2rh2(n2.theta), 1)))
    #                     log.system("\t\tDistance:", round(np.linalg.norm(n2.coords - n1.coords), 2))
    #                     count += 1
    #                 else:
    #                     log.system("WARNING: This crossover is between matching vectors. May indicate that the strand "
    #                           "directions of these adjacent modules is wrong.")
    #     log.system("Crossover count:", count)

    def validate_xover_distance(self):
        """
        Look for xover bonds that are far apart which may indicate that they were created at the wrong angle
        :return: printout
        """
        count = 0
        dist_list = []
        for plane in self.planes:
            for module in plane.modules:
                for xover in module.objxovers:
                    n1 = xover.n1
                    n1n5 = n1.n5
                    n1n3 = n1.n3
                    n2 = xover.n2
                    n2n5 = n2.n5
                    n2n3 = n2.n3

                    n1n5pos = n1n5.center + n1n5.nucl_vec
                    n1n3pos = n1n3.center + n1n3.nucl_vec
                    n2n5pos = n2n5.center + n2n5.nucl_vec
                    n2n3pos = n2n3.center + n2n3.nucl_vec

                    dist1 = round(np.linalg.norm(n1n5pos - n2n3pos), 2)
                    dist2 = round(np.linalg.norm(n1n3pos - n2n5pos), 2)
                    dist_list.append(dist1)
                    dist_list.append(dist2)
                    # log.system(n1n5.numid, n2n3.numid, dist1, n1n3.numid, n2n5.numid, dist2)

        log.system("\tMin: {} \n\tMax: {} \n\tMean: {} \n\tMedian: {} \n\tStd. Dev.: {}".format(
            min(dist_list),
            max(dist_list),
            np.mean(dist_list),
            np.median(dist_list),
            np.std(dist_list)
        ))
        return True

if __name__ == "__main__":
    pass