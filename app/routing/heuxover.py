"""
Created by Daniel Fu (Reif Lab, Duke University) at 12/4/2020

Name        : heuxover
Project     : cadaxisdna
Description :
Interpreter : Python 3.7.4
"""

from .helper import strandnav, log, mymath
from . import crossover
from .exceptions import *
from app import config
import random
import os


class CrossoverSet:
    """
    A collection of all crossover positions across the current Origami
    """

    def __init__(self, dnaorigami, **kwargs):
        """
        Initiates upon every edge
        """
        self.dnaorigami = dnaorigami
        if "exist_pos" in kwargs:
            exist_pos = kwargs.get("exist_pos")
            self.edges = list(exist_pos.keys())
            self.positions = {}
            for edge in self.edges:
                self.positions[edge] = [x for x in exist_pos[edge]]
        else:
            self.edges = dnaorigami.get_edges()
            self.positions = {}
            for edge in self.edges:
                self.positions[edge] = []

        # For internal tracking for producing crossover labels
        self.active_xovers = []

    def init_basic(self, bondlen=config.VALIDXOVERTHRESHBP,
                   xoverspacing1=config.VALIDXOVERSPACING_SAME,
                   xoverspacing2=config.VALIDXOVERSPACING_ADJ):
        """
        Sets a naive collection of proposed xover positions
        :param bondlen: see VALIDXOVERTHRESHBP
        :param xoverspacing1: valid xover spacing for same edge adjacent crossovers
        :param xoverspacing2: valid xover spacing for adjacent edge adjacent crossovers
        :return: None, sets self.positions
        """
        strand = 'stap'
        for edge in self.edges:
            module1 = edge.directed['from']
            module2 = edge.directed['to']
            xovers1 = module1.getxoversto(module2, 'stap')
            xovers2 = module2.getxoversto(module1, 'stap')

            # Filter by bond length
            proposed_xovers = crossover.filtervalid(xovers1, xovers2, thresh=bondlen)
            # Filter positions overlapping with existing xovers
            proposed_xovers = crossover.filteroverlap(proposed_xovers)
            # Filter positions too close to existing xovers
            proposed_xovers = self.filterdist2xover(proposed_xovers, xoverspacing1, xoverspacing2)
            # Filter duplicates in the list
            proposed_xovers = crossover.filterduplicate(proposed_xovers)
            # Ignores XOVERMODE, always dynamic, only ignores number of request xovers
            # Config, strand type, and the list of valid xovers must be non-empty
            # Ignore if list is empty (OK to initialize empty edge)
            if config.XOVERMODE == 'dynamic' and strand == 'stap' and proposed_xovers:
                proposed_xovers = crossover.sort_xovers(proposed_xovers)
                proposed_xovers = crossover.filterspacing(proposed_xovers)

            # Error flagging
            if not proposed_xovers:  # Empty
                log.log_warning("The crossover set between \n\n\t{}#{}\n\n\tand\n\n\t{}#{}\n\tis empty."
                                .format(str(module1), module1.up().up().get_module_index(module1),
                                        str(module2), module2.up().up().get_module_index(module2)))
            self.positions[edge] = proposed_xovers

    def filterdist2xover(self, xoverpairs,
                         same_edge_spacing=config.VALIDXOVERSPACING_SAME,
                         adj_edge_spacing=config.VALIDXOVERSPACING_ADJ,
                         same_edge=False):
        """
        Iterates through a list of XoverPair objects
        Filters for and then sorts by pairs that are farthest from any existing crossover
        Same as filterdist2feat but also works on unapplied XoverPairs
        :param xoverpairs: list of XoverPair
        :param same_edge_spacing: acceptable distance between adjacent crossovers spanning the same edge
        :param adj_edge_spacing: acceptable distance between adjacent crossovers spanning adjacent edges
        :param same_edge: whether to ignore adjacent edges
        :return: reduced list of XoverPairs
        """
        valid = []
        for xp in xoverpairs:
            same_edge_dist = self.distfromxover(xp, same_edge=True)
            if not same_edge:
                adj_edge_dist = self.distfromxover(xp, same_edge=False)
                if same_edge_dist >= same_edge_spacing and adj_edge_dist >= adj_edge_spacing:
                    valid.append(xp)
            else:
                if same_edge_dist >= same_edge_spacing:
                    valid.append(xp)
        return valid

    def refresh(self, bondlen=config.VALIDXOVERTHRESHBP,
                xoverspacing1=config.VALIDXOVERSPACING_SAME,
                xoverspacing2=config.VALIDXOVERSPACING_ADJ,
                scorethresh=config.HEUXOVER_RELAX_THRESH):
        """
        Do a single pass over all edges in a random order, selecting a new random collection of XoverPairs, then return
            self as a new CrossoverSet instance
        :param bondlen: see VALIDXOVERTHRESHBP
        :param xoverspacing1: valid xover spacing for same edge adjacent crossovers
        :param xoverspacing2: valid xover spacing for adjacent edge adjacent crossovers
        :param scorethresh: prioritizes the heuristic for largest regions
        :return: None, sets self.positions
        """
        # Iterate over the edges in a random order
        # Find all "regions" spanning adjacent helices on the same edge
        # For all long regions
        #   - Filter for valid crossovers
        #   - Filter for crossover distance on same edge
        #   - Randomly choose one crossover
        #   - Undo a crossover on any adjacent edge if it collides
        # print("Refreshing on {} {} {} {}".format(bondlen, xoverspacing1, xoverspacing2, scorethresh))
        edge_order = [e for e in self.positions]
        random.shuffle(edge_order)
        for edge in edge_order:
            module1 = edge.directed['from']
            module2 = edge.directed['to']
            if config.DYNAMIC_SPACING:  # STILL UNDER TEST, IGNORE
                min_bp = mymath.smaller_module(module1, module2).bp
                xoverspacing1 = min_bp / 21
                scorethresh = min_bp / 7
            try:  # Start at first XoverPair in list
                init_xp = self.positions[edge][0]
            except IndexError:  # If the edge is empty, add the full length of the module as a Xover insertion region
                if module1.bp >= module2.bp:
                    self.ins_xover(module1.helix.stap, edge, bondlen, xoverspacing1, xoverspacing2)
                else:
                    self.ins_xover(module2.helix.stap, edge, bondlen, xoverspacing1, xoverspacing2)
                # log.system("WARNING: XoverEdge list was empty, so no XoverPair was added.")
                continue
            next_xp = init_xp
            while True:  # Otherwise count from the first XoverPair around the module
                (region, next_xp) = self.region_to_next_xoverpair(next_xp)
                if len(region) >= scorethresh:
                    self.ins_xover(region, edge, bondlen, xoverspacing1, xoverspacing2)
                # Need to catch an error if init_xp is replaced during the ins_xover step
                if next_xp not in self.positions[self.get_xover_edge(next_xp)]:
                    # This should ONLY happen if the replacement is sharing nucleotides
                    # so find the XoverPair that shares nucleotides with next_xp
                    for xp in self.positions[self.get_xover_edge(next_xp)]:
                        if any(n in xp for n in next_xp.get_nucls()):
                            next_xp = xp
                            break
                    else:
                        raise RuntimeError("Replacement XoverPair does not intersect previous XoverPair.")

                # If fully traversed the module, exit the while loop
                if any(n in init_xp for n in next_xp.get_nucls()):
                    break

    def ins_xover(self, region, ref_edge,
                  bondlen=config.VALIDXOVERTHRESHBP,
                  xoverspacing1=config.VALIDXOVERSPACING_SAME,
                  xoverspacing2=config.VALIDXOVERSPACING_ADJ):
        """
        Inserts a crossover into an empty region
        :param region: list of Nucleotides to insert new crossover
        :param ref_edge: reference Edge for getting the 'to' and 'from' modules
        :param bondlen: crossover bond length
        :param xoverspacing1: spacing to same edge crossovers
        :param xoverspacing2: spacing to adjacent edge crossovers
        :return: inserts element in self.positions
        """
        module1 = ref_edge.directed['from']
        module2 = ref_edge.directed['to']
        lim_spacing_offset = 0
        xoverspacing1 += config.LIMSPACINGOFFSET
        xovers1 = module1.getxoversto(module2, 'stap')
        xovers2 = module2.getxoversto(module1, 'stap')
        bondlen = bondlen + (ref_edge.thresh - config.INTERHELICAL)

        # RUN LOOPS TO FIND A SET OF VALID CROSSOVERS
        while True:
            # Filter for valid crossovers
            xplist = crossover.filtervalid(xovers1, xovers2, thresh=bondlen)
            # Only keep those that completely fall in the region
            xplist = crossover.filtertoregion(region, xplist)
            if not xplist:
                log.system("WARNING: Increasing bond length locally from {1.1f} to {1.1f}".format(bondlen, bondlen+0.1))
                bondlen = round(bondlen + 0.1, 1)
                continue
            else:
                break
        while True:
            # Filter for crossover distance on same edge
            xplist = self.filterdist2xover(xplist, xoverspacing1, xoverspacing2, same_edge=True)
            if not xplist:
                lim_spacing_offset += 1
                if lim_spacing_offset > config.LIMSPACINGOFFSET:
                    # Give up if we've reached a limit
                    # Do not violate minimum spacing
                    break
                xoverspacing1 -= lim_spacing_offset
            else:
                break
        try:
            # print("Inserting crossover from {} choices into gap {} on {}".format(len(xplist), len(region), ref_edge))
            choose_xp = random.choice(xplist)
            # Must run clear_adj before add_xover because clear_adj can also look for overlaps
            # Otherwise it would immediately remove the added XoverPair
            self.clear_adj(choose_xp, int(xoverspacing2))
            self.add_xover(choose_xp)
        except IndexError:  # random.choice does not work for empty xplist
            # log.system("WARNING: Failed to insert crossover because no positions were available.")
            pass

    def clear_adj(self, x, r):
        """
        Removes other XoverPairs along any adjacent edges within the specified range
        :param x: source XoverPair
        :param r: int base pairs around init_xp
        :return: None, modifies self.positions
        """
        check_nucl = [x.n1.n5, x.n1.n3, x.n2.n5, x.n2.n3]
        for i in range(r):
            for n in check_nucl:
                n_xp = self.nucl_top_xover(n)
                if n_xp:
                    self.remove_xover(n_xp)
                    # print("Removed {} at range {} from {}".format(n_xp, i, x))
            check_nucl = strandnav.splash(*check_nucl)

    def distfromxover(self, x, same_edge=False, search_scaf=True):
        """
        Calculates distance to nearest XoverPair. Can also accomodate scaffold crossovers.
        :param x: XoverPair
        :param same_edge: search on same edge (True) or include adjacent edges (False)
        :param search_scaf: include scaffold strand XoverPairs in search?
        :return: int distance
        """
        init_nucl = [x.n1.n5, x.n1.n3, x.n2.n5, x.n2.n3]
        check_nucl = strandnav.splash(*init_nucl)
        check_nuclc = [n.Comp for n in check_nucl]
        d = 1
        while True:
            if any(n in init_nucl for n in check_nucl):
                if d == 1:
                    raise RuntimeError("distfromxover may have exited prematurely.")
                return d
            for n in check_nucl:  # Default: Check on staple strand
                nucl_xover = self.nucl_top_xover(n, detect_scaf=False)
                if nucl_xover:  # If it found an XoverPair
                    if strandnav.xover_common_edge(x, nucl_xover, same_edge):
                        return d
            if search_scaf:  # Check on scaffold strand
                for n in check_nuclc:
                    if strandnav.isxover(n):
                        if same_edge:
                            nucl_xover = self.dnaorigami.findxover(n)
                            if strandnav.xover_common_edge(x, nucl_xover, same_edge):
                                return d
                        else:
                            return d
            # Step outwards from source by 1 nucleotide each direction
            check_nucl = strandnav.splash(*check_nucl)
            check_nuclc = [n.Comp for n in check_nucl]
            d += 1

    def add_xover(self, item):
        """
        Adds the XoverPair into the collection
        :return: None
        """
        m1 = item.n1.n5.get_top_module()
        m2 = item.n2.n3.get_top_module()
        for edge in self.edges:
            if m1 in edge and m2 in edge:
                self.positions[edge].append(item)
                # log.system("Inserted {} on edge {}.".format(item, edge))
                return True
        raise RuntimeError("Addition of {} failed.".format(item))

    def remove_xover(self, item):
        """
        Removes the XoverPair from the collection
        :return: None
        """
        m1 = item.n1.n5.get_top_module()
        m2 = item.n2.n3.get_top_module()
        for edge in self.edges:
            if m1 in edge and m2 in edge:
                try:
                    self.positions[edge].remove(item)
                    # log.system("Removed {} from edge {}.".format(item, edge))
                except ValueError:
                    log.system("[WARNING] {} not among proposed crossover positions, "
                                         "possibly already removed.".format(item))

    def region_to_next_xoverpair(self, init_xp):
        """
        Returns the region from an XoverPair to the next XoverPair
        :param init_xp: XoverPair
        :return: list of Nucleotides, next XoverPair
        """
        # Initial and end nucl for traversal
        e1 = init_xp.n1.n5
        e2 = init_xp.n2.n3

        # Traverse on edge 1 first, then edge 2, then use the longer strand
        # Traverse edge 1
        e1 = e1.__strand5__
        e1_nucl = [e1]
        while True:  # Traverse edge 1 and track Nucleotides until reaching the next XoverPair
            e1_xp = self.nucl_top_xover(e1)  # Check for a corresponding XoverPair
            if e1_xp:  # Found a corresponding XoverPair
                if strandnav.xover_common_edge(init_xp, e1_xp, True):  # Has to span the same modules
                    break
            e1 = e1.__strand5__  # Continue traversing
            e1_nucl.append(e1)
        # Traverse edge 2
        e2 = e2.__strand3__
        e2_nucl = [e2]
        while True:  # Traverse edge 2 and track Nucleotides until reaching the next XoverPair
            e2_xp = self.nucl_top_xover(e2)  # Check for a corresponding XoverPair
            if e2_xp:  # Found a corresponding XoverPair
                if strandnav.xover_common_edge(init_xp, e2_xp, True):  # Has to span the same modules
                    break
            e2 = e2.__strand3__  # Continue traversing
            e2_nucl.append(e2)
        if e2_xp != e1_xp:
            log.system("XoverPair1: {}\n XoverPair2: {}\n"
                       "Started from {}\n\n"
                       "Regions:\n"
                       "e1 5' {}\n"
                       "e2 3' {}\n".format(e1_xp, e2_xp, init_xp, [n.numid for n in e1_nucl],
                                           [n.numid for n in e2_nucl]))
            # log.system("All crossover on edges (alt. method): {}".format(self.positions[self.get_xover_edge(init_xp)]))
            # Error check to make sure it is the same xover detected in the previous loop
            raise RuntimeError("Traversing module 1 and module 2 did not encounter the same crossover.")

        # Use the longer strand
        if len(e1_nucl) >= len(e2_nucl):
            return e1_nucl, e1_xp
        else:
            return e2_nucl, e2_xp

    def score(self):
        """
        Calculate and return the crossover spacing score of the proposed xovers
        :return: int score
        """
        crossover_scores = []
        for edge in self.positions:
            try:  # Start at first XoverPair in list
                init_xp = self.positions[edge][0]
            except IndexError:  # If the edge is empty, add the full length of the module as a gap
                module1 = edge.directed['from']
                module2 = edge.directed['to']
                crossover_scores.append(max([module1.bp, module2.bp]))
                continue
            next_xp = init_xp
            while True:  # Otherwise count from the first XoverPair around the module
                (region, next_xp) = self.region_to_next_xoverpair(next_xp)
                # Add the gap score
                crossover_scores.append(len(region) + 1)
                # If fully traversed the module, exit the while loop
                if next_xp == init_xp:
                    break
        return max(crossover_scores)

    def count(self):
        """
        Returns the number of currently held potential crossover positions
        :return: int count
        """
        count = 0
        for edge in self.positions:
            for xp in self.positions[edge]:
                count += 1
        return count

    def print_xoverlabels(self, filename):
        """
        Integrates with heuxover.optimize_crossovers to save intermediate steps of the heuristic
        :param filename: string, file prefix to save as
        :return: None
        """
        xovers_filename = os.path.join(self.dnaorigami.outputdir, filename + '.xovers')
        os.makedirs(os.path.dirname(xovers_filename), exist_ok=True)

        # Apply, save, then undo
        self.apply()
        self.dnaorigami.print_oxDNA(filename=filename, print_crossovers=True)
        self.undo()

    def apply(self):
        """
        Applies all XoverPairs found in self.position into FullXovers
        :return:
        """
        for edge in self.positions:
            for xp in self.positions[edge]:
                fx = xp.apply()
                self.dnaorigami.protect(xp.n1.to_tuple[0], config.PROTECT)
                self.dnaorigami.protect(xp.n2.to_tuple[0], config.PROTECT)
                fx.settype('stap')
                xp.n1.up().addxover(fx)
                xp.n2.up().addxover(fx)
                self.active_xovers.append(fx)

    def undo(self):
        if not self.active_xovers:
            raise RuntimeError("Crossover set cannot be reversed as there are no activated crossovers.")
        else:
            for fx in self.active_xovers:
                fx.undo()
        self.active_xovers = []


    def nucl_top_xover(self, nucl, detect_scaf=False):
        """
        Checks whether the Nucleotide coordinates to an XoverPair
        :param nucl: Nucleotide
        :param detect_scaf: Whether to run the scaffold snippet
        :return: Return the XoverPair, else return False
        """
        for edge in self.positions:
            for x1 in self.positions[edge]:
                if nucl in x1:
                    return x1
        if detect_scaf:
            # Test for scaffold crossover, but if found, still return a virtual staple strand XoverPair as reference
            try:  # At this point, there should only be scaffold xovers
                xp = self.dnaorigami.findxover(nucl.Comp).complement()
                return xp
            except NuclNotFound:
                pass
        return False

    def get_xover_edge(self, xp):
        """
        Returns the edge that the XoverPair spans
        :param xp: XoverPair
        :return:
            - path.Edge if found
            - bool False otherwise
        """
        m1 = xp.n1.n5.get_top_module()
        m2 = xp.n2.n3.get_top_module()
        for edge in self.edges:
            if m1 in edge and m2 in edge:
                if xp.spans(m1, m2):
                    return edge
        return False


class Mixin:
    def optimize_crossovers(self, save_steps=False):
        """
        Controller method for doing simulated annealing optimization of crossovers
        :return: None
        """
        # temperature = round(max([m.bp for m in self.get_modules()])/5)
        temperature = config.SA_INIT_TEMP
        fails = 0
        step = 0
        while True:
            current_score = self.xover_set.score()
            # Generate a new instance of crossovers
            new_step = CrossoverSet(self, exist_pos=self.xover_set.positions)
            new_step.refresh(bondlen=config.VALIDXOVERTHRESHBP,
                             xoverspacing1=config.VALIDXOVERSPACING_SAME,
                             xoverspacing2=config.VALIDXOVERSPACING_ADJ,
                             scorethresh=config.HEUXOVER_RELAX_THRESH)
            # Check score
            new_score = new_step.score()
            # Prints the current score
            log.system("Testing new config of score {}".format(new_score))
            # Prints the number of crossovers
            log.system("Current crossover count: {}".format(new_step.count() + self.num_xovers()))
            # Saves the current configuration
            if save_steps:
                log.system("Exporting step")
                new_step.print_xoverlabels("step"+"{:03}".format(step))
                step += 1
            if new_score <= current_score + temperature:
                log.system("Moving from {} to {} at temperature {}".format(current_score, new_score, temperature))
                fails = 0
                self.xover_set = new_step
                if new_score <= config.SA_ACCEPT:
                    break
                temperature = temperature - config.SA_RATE
            elif temperature < 3.0:  # Taper down on consecutive failures if we are nearing minimum
                temperature -= config.SA_RATE/10
            else:
                if fails > config.SA_FAIL_COUNT:
                    temperature -= config.SA_RATE
                else:
                    fails += 1
            if temperature < 0:
                break
        self.xover_set.apply()


if __name__ == "__main__":
    pass
