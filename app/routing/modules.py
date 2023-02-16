# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 21:23:56 2018

@author: Dan

CHANGELOG:
    7/12/2018:
        - Write RingModule class
    8/20/2018:
        - Write Plane class
    9/12/2018:
        - Renamed InitHelix to DNAHelix
    9/29/2018:
        Added new attributes to RingModule to help track crossovers
            methods: addxover, delxover
            values: objxovers (so as to not coincide with xovers)
"""
from .helper import mytrig, mymath, rmatrix, dnaconnector, log, strandnav
from scipy.spatial.transform import Rotation as R
from .score import seeding
from . import dnahelix, crossover, nicking, calibrate, makeshape
from .exceptions import DistanceError
from app import config
import numpy as np


def fitarc(xf, bp):
    """
    Calculates a tbx and numxovers to fit the input xover factor and bp circumference
    :param xf:
    :param bp:
    :return:
    """
    tbx = 2
    numxovers = xf

    bpt = bp / numxovers / tbx
    if bpt < 9:
        log.system("ERROR: Extremely small ring possibly detected. Cannot fit lower bound of crossovers.")
        raise ValueError

    while bpt < 9 or bpt > 12:
        tbx += 1
        if tbx > 5:
            tbx = 2
            numxovers += xf
        bpt = bp / numxovers / tbx
    return tbx, numxovers

class Plane(object):
    """
    Plane is a z-level cross-section of all modules at that height
    """

    def __init__(self, origin, pitch=0, yaw=0):
        """
        :param center: reference origin for the plane
        :param xyz: surface normal
        """
        # Default Axes
        self.bx = np.array([1, 0, 0])
        self.by = np.array([0, 1, 0])
        self.bz = np.array([0, 0, 1])

        r = R.from_euler('xyz', [yaw, 0, pitch], degrees=True)
        rmat = r.as_dcm()

        self.bx = np.matmul(rmat, self.bx)
        self.by = np.matmul(rmat, self.by)
        self.bz = np.matmul(rmat, self.bz)
        self.basis = [self.bx, self.by, self.bz]
        log.debug(__name__, "Creating PLANE (h={})".format(origin[2]))
        self.modules = []
        self.origin = origin
        self.height = origin[2]
        self.adjacent = []
        self.yaw = yaw
        self.pitch = pitch
        self.top = None

    def add_circle(self, idx, bp, dir_bit, **kwargs):
        """
        Adda a circle onto the plane
        :param idx: current nucleotide index
        :param bp: base pair circumference
        :param dir_bit: 5'-3' direction of the helix
        :return:
        """
        if config.TWIST_NORMALIZED:
            new_geom = mymath.min_to_bpt(bp)
            turns = new_geom['turns']
            bp += new_geom['i']
        else:
            try:
                turns = kwargs['tbx'] * kwargs['numxovers']
            except KeyError:
                log.system("WARNING: Input was not given turns per crossovers and number of crossovers for calculating "
                           "turns. Refitting for those parameters based on XOVER_FACTOR.")
                (tbx, numxovers) = fitarc(config.XOVER_FACTOR, bp)
                turns = tbx * numxovers
                kwargs['tbx'] = tbx
                kwargs['numxovers'] = numxovers

        shape_obj = makeshape.CircleOnPlane(bp, self.origin, self.basis)
        dna_obj = RingModule(shape_obj.num_points, turns, dir_bit, self.origin, self.bz, **kwargs)
        dna_obj.applyshape(shape_obj, idx)
        dna_obj.top = self  # Set hierarchy relationship
        self.modules.append(dna_obj)
        idx += int(dna_obj.shape.num_points * 2)
        return dna_obj, idx

    def add_arc(self, idx, numpts, dir_bit, a1, arclen, center):
        """
        Adds an arc onto the plane
        :param idx: current nucleotide index
        :param numpts: base pair length
        :param dir_bit: 5'-3' direction of the helix
        Use three-point arc drawing:
        :param a1: first poitn on arc
        :param arclen: arc length
        :param center: center of arc
        :param forced_circumference: Fixes the length of the arc
        :param align: Makes the nucleotides at the arc's endpoints have a specified helical angle
        :return:
        """
        new_geom = mymath.min_to_bpt(numpts)
        turns = new_geom['turns']
        numpts += new_geom['i']
        shape_obj = makeshape.ArcOnPlane(center, numpts, a1, arclen, self.origin, self.basis)
        dna_obj = ArcModule(shape_obj.num_points, turns, dir_bit, self.origin, self.bz)
        dna_obj.applyshape(shape_obj, idx)
        dna_obj.top = self  # Set hierarchy relationship
        self.modules.append(dna_obj)
        idx += int(dna_obj.shape.num_points * 2)
        return dna_obj, idx

    def add_line(self, idx, start, duvec, numpts, dir_bit, **kwargs):
        offset = kwargs.get("offset", 0)
        turns = numpts / 10.5
        shape_obj = makeshape.LineOnPlane(numpts, start, duvec)
        dna_obj = LineModule(shape_obj.num_points, turns, dir_bit, self.origin, self.bz)
        dna_obj.applyshape(shape_obj, idx, offset=offset)
        dna_obj.top = self  # Set hierarchy relationship
        self.modules.append(dna_obj)
        idx += int(dna_obj.shape.num_points * 2)
        return dna_obj, idx

    def add_module(self, obj):
        """
        Adds a module to the plane
        :param obj: module object
        :return: None
        """
        if all(z.coordinates[2] != self.height for z in obj.shape.graph):
            # Error handler if object is not in plane
            log.system([z.coordinates[2] for z in obj.shape.graph])
            raise ValueError("Cannot add new module because it is at height {} "
                             "but plane is at height {}".format(obj.height, self.height))
        else:
            self.modules.append(obj)
            obj.top = self  # Set hierarchy relation
            # self.top.last_idx += int(obj.shape.num_points * 2)  # Update idx

    def set_height(self, h):
        for module in self.modules:
            if module.height != h:
                log.out(__name__,
                        "Trying to set height {}\
                        but module r:{},bp:{}\
                        is at height {}".format(self.height, module.radius, module.bp, module.height))
                s = None
                while s != "Y" or s != "N":
                    s = input("Confirm change? Only use (Y/N)")
                if s == "Y":
                    module.height = h
                elif s == "N":
                    pass
        self.height = h

    def contains(self, obj):
        if obj in self.modules:
            return True
        return False

    def settop(self, top):
        """
        Sets the top level object in the DNA structure hierarchy
        :param top: object containing this object
        :return: None
        """
        self.top = top

    def up(self):
        """
        Returns the top level object in the DNA structure hierarchy
        :return: object containing this object, Origami(object)
        """
        return self.top


""" DNA helix shape modules """


def connect(module1, module2, angle=None, match=True):
    """
    connects two modules together at an automatically found endpoint.
    a valid endpoint is either the 5' and 3' ends of either module, respectively
    no valid endpoint is evaluated if both pairs of points are too far from each other
    :param module1: 1st module
    :param module2: 2nd module
    :param angle: offset angle by which these modules should be connected
    :param match: whether to match twist of adjcent modules or not
    :return: None
    """
    # Due to invariant, first index of sequence list is 5' end of strand and last index is 3'
    m1n3 = module1.helix.scaf[-1]
    m2n5 = module2.helix.scaf[0]
    # Error handler
    # Calculate the distance between pairs
    d = np.linalg.norm(m1n3.center - m2n5.center)
    if d > 1.0:
        log.log_warning("Module connections exited early. Can't connect these two modules. "
                        "They may not be adjacent or have opposing 5'-3' "
                        "directions. The separation distance of their 5' and 3' end is {}.".format(d))
        # debug.highlight_nucl(m1n3, mult=10)
        # debug.highlight_nucl(m2n5, mult=10)
        raise DistanceError
    # Calibrate their connection angle (optional) and connect at manually input or auto angle
    if match:
        if angle:
            calibrate.matchtwist(m1n3, m2n5, angle=angle)
        else:
            calibrate.matchtwist(m1n3, m2n5)
    dnaconnector.nuclconnect(m1n3, m2n5)
    dnaconnector.nuclconnect(m2n5.Comp, m1n3.Comp)
    dnaconnector.trueconnect(m1n3, m2n5)
    dnaconnector.trueconnect(m2n5.Comp, m1n3.Comp)


def relax_junctions(refnucl):
    """
    Not yet implemented.
    Propogates an angular difference between two adjacent nucleotides into adjacent pairs along the strand to
        try and even out the twist
    :param modules:
    :return:
    """
    strand = strandnav.getstrand(refnucl)
    for nucl in strand:
        nucl3 = nucl.toThree
        try:
            t3 = mytrig.acute_diff(nucl3.theta, nucl.theta)
        except AttributeError:
            pass
        if t3 > 40:
            pass


class JoinedStrand:
    """
    DNA double helix module of any shape created by combining multiple independent modules into a single super-module
    Combines strands of piecewise modules into a single strand
    """
    def __init__(self, piecewise_modules):
        self.modules = []
        self.adjacent = []
        if not all(m.note_helixnum == piecewise_modules[0].note_helixnum for m in piecewise_modules):
            log.log_error("Modules may not all be part of the same strand.")
        for m in piecewise_modules:
            m.topjs = self
            for scaffold_nucl in m.helix.scaf:
                scaffold_nucl.strand = self
            for staple_nucl in m.helix.stap:
                staple_nucl.strand = self
            self.modules.append(m)

    def set_adjacent(self, strand):
        if strand not in self.adjacent:
            self.adjacent.append(strand)

    def get_length(self):
        total_length = 0
        for m in self.modules:
            total_length += m.bp
        return total_length

    def get_staple(self):
        return strandnav.absgetstrand(self.modules[0].helix.stap[0])

    def get_scaffold(self):
        return strandnav.absgetstrand(self.modules[0].helix.scaf[0])

    def crossovers_to(self, js):
        """
        Returns the FullXover objects pointed to other JoinedStrands in order of sequential appearance
        :param js: target joined strand
        :return: List of Xovers, in order of appearance
        """
        if js not in self.adjacent:
            return []

        # Get a reference nucleotide for the entire strand
        ref_nucl = self.modules[0].helix.stap[0]
        this_strand = strandnav.absgetstrand(ref_nucl)
        last_xover = None
        init_idx = 0
        # Cycle to a crossover
        for idx in range(len(this_strand)):
            curr_nucl = this_strand[idx]
            next_nucl = this_strand[(idx + 1) % len(this_strand)]
            # Query the crossovers
            try:
                x1 = curr_nucl.get_top_module().get_xover(curr_nucl)
            except ValueError:
                try:
                    x1 = curr_nucl.get_top_module().get_xover(curr_nucl.Comp)
                except ValueError:
                    x1 = None
            try:
                x2 = next_nucl.get_top_module().get_xover(next_nucl)
            except ValueError:
                try:
                    x2 = next_nucl.get_top_module().get_xover(next_nucl.Comp)
                except ValueError:
                    x2 = None
            # Crossovers exist, are the same crossover and spans the input strands
            if (x1 and x2) and (x1 is x2) and (x1.spans(self, js) and x2.spans(self, js)):
                init_idx = (idx + 1) % len(this_strand)
                last_xover = x1
                break

        # Cycle the list of crossovers to begin at the latter half of it so it is not counted twice
        helix_xovers = mymath.cycle_list(this_strand, init_idx)

        crossovers = []
        # Go through every base in the strand
        for base in helix_xovers:
            # Only check crossovers
            if strandnav.isxover(base) or strandnav.isxover(base.Comp):
                # checks for either staple or scaffold crossover
                try:
                    xobj = base.get_top_module().get_xover(base)
                except ValueError:
                    xobj = base.get_top_module().get_xover(base.Comp)
                if xobj == last_xover:  # Ignore consecutive bases that are part of same crossover
                    pass
                elif xobj.spans(self, js):  # Spans the input strands, not yet counted
                    if xobj not in crossovers:
                        crossovers.append(xobj)
                    last_xover = xobj  # Reset last reference point
                else:
                    pass
        return crossovers

    def __repr__(self):
        return self.modules[0].note


class Shape:
    """
    A list of nucleotides following a given shape trace
    """

    class Note:
        def __init__(self):
            pass

        def __get__(self, instance, owner):
            return instance.note_string

        def __set__(self, instance, value):
            instance.note_string = value
            s = mymath.delete_chars(value, [",", "(", ")"])
            s = s.split()
            instance.note_helixnum = (int(s[0]), int(s[1]))
            instance.note_section = str(s[2])

    note = Note()

    def __init__(self, bp, turns, dirBit, origin, normal):
        # calculate and set properties
        self.bp = bp  # base pair circumference
        self.turns = turns  # number of full turns counted across the entire helix
        self.height = origin[2]
        self.origin = origin
        self.normal = normal
        self.dirBit = dirBit
        self.bpt = self.bp / self.turns  # bases per turn
        self.apb = 360.0 / self.bpt  # angle per base

        # the adjacent list must go to the same shape that is positioned in parallel to the current shape
        self.adjacent = {}  # adjacency list

        self.shape = None
        self.helix = None
        # TODO: Consider making this a set
        self.objxovers = []
        self.top = None
        self.topjs = None

        # gives the object an ID number
        self.numid = -1

        # Adds a label for debugging purposes
        self.note_string = "(-1, -1) X"
        self.note_helixnum = (-1, -1)
        self.note_section = "X"

    def report_breaks(self):
        """
        DEBUG method
        :return:
        """
        for nucl in self.helix.scaf:
            if nucl.__strand3__ == -1:
                log.system("DEBUG: Dangling 3' scaf strand found: {} {}".format(nucl.numid, nucl.get_top_module()))
            elif nucl.__strand5__ == -1:
                log.system("DEBUG: Dangling 5' scaf strand found: {} {}".format(nucl.numid, nucl.get_top_module()))

        for nucl in self.helix.stap:
            if nucl.__strand3__ == -1:
                log.system("DEBUG: Dangling 3' stap strand found: {} {}".format(nucl.numid, nucl.get_top_module()))
            elif nucl.__strand5__ == -1:
                log.system("DEBUG: Dangling 5' stap strand found: {} {}".format(nucl.numid, nucl.get_top_module()))

    def clean(self):
        """
        Looks into its helix for short helices that are unconnected to the main structure
        :return: None, removes excess nucleotides from Origami
        """
        remove_buffer = []
        for nucl in self.helix.scaf:
            if nucl.__strand5__ == -1:
                this_nucl = nucl
                while this_nucl.__strand3__ != -1:
                    remove_buffer.append(this_nucl)
                    this_nucl = this_nucl.__strand3__
                    if len(remove_buffer) > self.bp:
                        raise RuntimeError("Strand has looped without finding the other noncyclic terminus.")
                remove_buffer.append(this_nucl) # Get the last position
        if remove_buffer:
            # print("DEBUG: Length before {}".format(len(self.helix.scaf)))
            for nucl in remove_buffer:
                self.helix.scaf.remove(nucl)
            # print("DEBUG: Removed {} from scaffold on {}.".format(len(remove_buffer), self))
            # print("DEBUG: Length after {}".format(len(self.helix.scaf)))

        remove_buffer = []
        for nucl in self.helix.stap:
            if nucl.__strand5__ == -1:
                this_nucl = nucl
                while this_nucl.__strand3__ != -1:
                    remove_buffer.append(this_nucl)
                    this_nucl = this_nucl.__strand3__
                    if len(remove_buffer) > self.bp:
                        raise RuntimeError("Strand has looped without finding the other noncyclic terminus.")
                remove_buffer.append(this_nucl)  # Get the last position
        if remove_buffer:
            # print("DEBUG: Length before {}".format(len(self.helix.stap)))
            for nucl in remove_buffer:
                self.helix.stap.remove(nucl)
            # print("DEBUG: Removed {} from staple on {}.".format(len(remove_buffer), self))
            # print("DEBUG: Length after {}".format(len(self.helix.stap)))

    # Topology methods
    def applyshape(self, shapeobj, init_base_index, align=False, **kwargs):
        """
        Applies nucleotides onto the list of points describing a shape
        :param shapeobj:
        :param init_base_index:
        :return:
        """
        self.shape = shapeobj

        if not align:
            init_angle = 90 + self.dirBit * 240
        else:
            # self.apb = self.apb + 2 * self.apb/self.bp
            init_angle = (align - self.apb/2) % 360

        if 'offset' in kwargs:
            init_angle = kwargs.get('offset', 0)

        self.helix = dnahelix.DNAHelix(
            self.shape,
            self.apb,
            init_angle,
            init_base_index,
            self.dirBit, self)

    def settop(self, top):
        """
        Sets the top level object in the DNA structure hierarchy
        :param top: object containing this object
        :return: None
        """
        self.top = top

    def up(self):
        """
        Returns the top level object (Plane) in the DNA structure hierarchy
        :return: Plane
        """
        return self.top

    def get_top_plane(self):
        """
        :return: Topology containing Plane
        """
        return self.up()

    def get_top_origami(self):
        """
        :return: Topology containing Origami
        """
        return self.get_top_plane().up()

    def get_top_strand(self):
        """
        :return: Topology containing JoinedStrand
        """
        return self.topjs

    # Crossover methods
    def delxover(self, fx):
        self.objxovers.remove(fx)

    def addxover(self, fx):
        self.objxovers.append(fx)

    def get_xover(self, nucl):
        """
        Returns an xover object by a reference nucleotide
        :param nucl: Nucleotide
        :return: FullXover
        """
        # Look on this module
        for x in self.objxovers:
            if x.has_nucl(nucl):
                return x

        # Some debugging stuff, remove later
        # numids = []
        # for x in self.objxovers:
        #     [numids.append(nid) for nid in list(x.n1.numid)+list(x.n2.numid)]
        # print("Nucleotides in objxovers:", numids)
        # print("Is it an xover,", strandnav.isxover(nucl))
        # print("Is comp an xover,", strandnav.isxover(nucl.Comp))
        # print("nucl.numid", nucl.numid)
        # print(nucl.get_top_module())
        # print("nucl.Comp.numid", nucl.Comp.numid)
        # print(nucl.Comp.get_top_module())
        # print("nucl.toThree.numid", nucl.toThree.numid)
        # print(nucl.toThree.get_top_module())
        # print("nucl.toFive.numid", nucl.toFive.numid)
        # print(nucl.toFive.get_top_module())
        # print("nucl.__strand3__.numid", nucl.__strand3__.numid)
        # print(nucl.__strand3__.get_top_module())
        # print("nucl.__strand5__.numid", nucl.__strand5__.numid)
        # print(nucl.__strand5__.get_top_module())
        # print("nucl.Comp.toThree.numid", nucl.Comp.toThree.numid)
        # print(nucl.Comp.toThree.get_top_module())
        # print("nucl.Comp.toFive.numid", nucl.Comp.toFive.numid)
        # print(nucl.Comp.toFive.get_top_module())
        # print("nucl.Comp.__strand3__.numid", nucl.Comp.__strand3__.numid)
        # print(nucl.Comp.__strand3__.get_top_module())
        # print("nucl.Comp.__strand5__.numid", nucl.Comp.__strand5__.numid)
        # print(nucl.Comp.__strand5__.get_top_module())
        # log.log_warning("Module {} has no crossover corresponding to nucleotide {}.".format(self, nucl.numid))

        nucl_module = nucl.get_top_module()
        nucl3_module = nucl.__strand3__.get_top_module()
        nucl5_module = nucl.__strand5__.get_top_module()
        for x in nucl5_module.objxovers:
            if x.has_nucl(nucl):
                return x
        for x in nucl3_module.objxovers:
            if x.has_nucl(nucl):
                return x
        for x in nucl_module.objxovers:
            if x.has_nucl(nucl):
                return x
        raise ValueError("Crossover cannot be found.")

        # There is an exceptional case where if a crossover straddles a module junction within the strand,
        #   the crossover may be on an adjacent module
        # Look for the crossover in the entire structure and flag a warning
        # Conditional: If there is xover at this nucl but it was not detected in the module's objxover list
        # if strandnav.isxover(nucl) or strandnav.isxover(nucl.Comp):
        #     origami = nucl.get_top_origami()
        #     log.log_warning("Had to look for a corresponding crossover to {} in the overall structure."
        #                     .format(nucl.numid))
        #     try:
        #         return origami.findxover(nucl)
        #     except NuclNotFound:
        #         return origami.findxover(nucl.Comp)

    def getxoversto(self, to_ring, choose_strand):
        """
        Iterates through all bases on the strand
        Always goes 5' to 3'
        Returns an ordered list of NuclPair tuples
        Formatted (5',3') base of the valid crossover position
        :param to_ring:
        :param choose_strand:
        :return:
        """
        if choose_strand == 'scaf':
            try:
                return self.adjacent[to_ring]['scxovers']
            except KeyError:
                raise KeyError("{} not found in scaffolds adjacency list of {}.".format(to_ring, self))
        elif choose_strand == 'stap':
            try:
                return self.adjacent[to_ring]['stxovers']
            except KeyError:
                raise KeyError("{} not found in staples adjacency list of {}.".format(to_ring, self))
        else:
            raise ValueError("No strand {}. Please choose 'scaf' or 'stap'.".format(choose_strand))

    def setadjacent(self, module, angle):
        """
        Save valid scaf and stap crossovers
        :param module: Target module
        :param angle: float Angle of connection
        :return: None
        """
        # init key with empty xovers list value
        self.adjacent[module] = {'scxovers': [], 'stxovers': []}
        self.setxoverstoangle(angle, module, 'stap')
        self.setxoverstoangle(angle, module, 'scaf')

        log.system("{} set adjacent to {}".format(module, self))

    def avg_angle_to(self, module):
        """
        Finds the angle to an adjacent module as an average of all its nearest neighbors per base
        :param module: modules.Shape
        :return: int angle value
        """
        strand = self.helix.getstrand('stap')
        all_angles = []
        for base in strand:
            all_angles.append(mymath.base_to_module_angle(base, module) % 360)

        return sum(all_angles)/len(all_angles)

    def shortest_dist_to(self, module):
        """
        Finds the shortest distance to an adjacent module
        :param module: modules.Shape
        :return: int distance value
        """
        strand = self.helix.getstrand('stap')
        all_dists = []
        for base in strand:
            all_dists.append(mymath.base_to_module_dist(base, module))

        return min(all_dists)

    def avg_dist_to(self, module):
        """
        Finds the distance to an adjacent module as an average of all its nearest neighbors per base
        :param module: modules.Shape
        :return: int distance value
        """
        strand = self.helix.getstrand('stap')
        all_dists = []
        for base in strand:
            all_dists.append(mymath.base_to_module_dist(base, module))

        return sum(all_dists) / len(all_dists)

    def setxoverstoangle(self, angle, module, choose_strand):
        """
        Finds pairs of adjacent nucleotides that may be valid crossover positions to the target module
        :param angle: The angle to the target module
            4/4/2020 Deprecated for spherical coordinates method
        :param module: Target module
        :param choose_strand: Chosen strand ('scaf' or 'stap')
        :return: None, calls self.setxoverslist to modify the self.adjacent attribute
        """
        xovers = []
        strand = self.helix.getstrand(choose_strand)
        # Go through each pair of adjacent nucleotides in the strand
        for base in strand:
            base1 = base
            base2 = base.__strand3__
            target_angle = mymath.base_to_module_angle(base, module)
            # target_vector = mymath.base_to_module_vector(base, module)
            try:
                # Get the theta of the Nucleotide's, already on the correct basis
                # Base2 may have a slightly different basis, if module is curved, but it shouldn't be so different that
                # it would introduce significant error, so it is not normalized to the base1 basis.
                base1theta = base1.theta
                base2theta = base2.theta
                flag_angle = mytrig.angle_between(base1theta, target_angle, base2theta)

                # base1vector = base1.vector
                # base2vector = base2.vector
                # bx = base1.bx
                # by = base1.by

                # Project everything to base1's plane
                # base2vector_p = mymath.project_3dvector_to_2dbasis(base2vector, bx, by)
                # target_vector_p = mymath.project_3dvector_to_2dbasis(target_vector, bx, by)
                #
                # flag_angle = mytrig.vector_between(base1vector, target_vector_p, base2vector_p)

                # DEBUGGING
                # print("Checking NuclPair({}, {}) for a crossover location 5': {} t: {} 3': {} | {}".format(
                #     base1.numid, base2.numid,
                #     round(base1theta, 2), round(target_angle, 2), round(base2theta, 2), flag_angle))
                if flag_angle:
                    xovers.append(crossover.NuclPair(base1, base2))
            # Error handler if base.toThree is a break, thus type 'int' and no attribute theta
            except AttributeError:
                pass
        # Catch error if no xover positions were found.
        if not xovers:
            raise RuntimeError("No valid crossover positions were set.")
        # Add them to respective dictionaries
        self.setxoverslist(xovers, module, choose_strand)

    def setxoverslist(self, xoverslist, module, choose_strand):
        if choose_strand == 'scaf':
            self.adjacent[module]['scxovers'] = xoverslist
        elif choose_strand == 'stap':
            self.adjacent[module]['stxovers'] = xoverslist
        else:
            raise ValueError("No strand {}. Please choose 'scaf' or 'stap'.".format(choose_strand))

    # EXPORT METHODS
    # Sets the internal variables for each nucleotide normalized to a caDNAno lattice
    def set_cadnano(self, size, num):
        self.cdna_num = num
        blank = size - self.bp
        padding = blank / 2
        if not padding.is_integer():
            raise ValueError("Padding is not an integer value.")

        # check if any break or crossover is on a border
        # conditions:
        while strandnav.isfeature(self.helix.scaf[0]) or \
                strandnav.isfeature(self.helix.scaf[-1]) or \
                strandnav.isfeature(self.helix.stap[0]) or \
                strandnav.isfeature(self.helix.stap[-1]):
            # action:
            self.helix.scaf = mymath.cycle_list(self.helix.scaf, 1)
            self.helix.stap = mymath.cycle_list(self.helix.stap, 1)

        # sets the internal nucleotide variables
        for ind, nucl in enumerate(self.helix.stap, start=0):
            nucl.cdna_num = num
            nucl.cdna_ind = int(ind + padding)
        for ind, nucl in enumerate(self.helix.scaf, start=0):
            nucl.cdna_num = num
            nucl.cdna_ind = int(ind + padding)

    # set each ring into the 4-tuple format of each caDNAno grid square
    def format_cadnano(self, size):
        all_stap_colors = [16204552, 11184640, 12060012, 13369344, 8947848, 243362, 16225054, 7536862, 3355443,
                           29184, 5749504]
        # [F74308,AAAA00,B8056C,CC0000,888888,03B6A2,F7931E,7300DE,333333,007200,57BB00]
        # Red Orange, Yellow Green, Magenta, Red, Grey, Teal, Orange, Purple, Dark Grey, Dark Green, Bright Green
        blank = size - self.bp
        padding = blank / 2
        ostp = [[-1, -1, -1, -1] for p in range(int(padding))]
        ostn = []
        stap_colors = []
        for ind, nucl in enumerate(self.helix.stap, start=0):
            try:
                to5num = nucl.toFive.cdna_num
                to5ind = nucl.toFive.cdna_ind
            except AttributeError:
                to5num = -1
                to5ind = -1
            try:
                to3num = nucl.toThree.cdna_num
                to3ind = nucl.toThree.cdna_ind
            except AttributeError:
                to3num = -1
                to3ind = -1
            if nucl.toFive == -1:
                # random colors
                #                random.shuffle(all_stap_colors)
                #                stap_colors.append([ind+padding,all_stap_colors[0]])

                # has_seed coloring
                # red for no seed
                # grey for yes seed
                if seeding.has_seed(strandnav.getstrand(nucl), 14):
                    stap_colors.append([ind + padding, 8947848])
                else:
                    stap_colors.append([ind + padding, 13369344])

            #                # short region coloring
            #                # red for risk region < 5 nt
            #                # grey for OK
            #                if not region.has_short_region(strandnav.getstrand(nucl),4):
            #                    stap_colors.append([ind+padding,8947848])
            #                else:
            #                    stap_colors.append([ind+padding,13369344])

            ostn.append([to5num, to5ind, to3num, to3ind])

        oscp = [[-1, -1, -1, -1] for p in range(int(padding))]
        oscn = []
        for nucl in self.helix.scaf:
            try:
                to5num = nucl.toFive.cdna_num
                to5ind = nucl.toFive.cdna_ind
            except AttributeError:
                to5num = -1
                to5ind = -1
            try:
                to3num = nucl.toThree.cdna_num
                to3ind = nucl.toThree.cdna_ind
            except AttributeError:
                to3num = -1
                to3ind = -1
            oscn.append([to5num, to5ind, to3num, to3ind])

        ost = ostp + ostn + ostp
        osc = oscp + oscn + ostp

        return {'scaf': osc, 'stap': ost, 'colors': stap_colors}

    # QUERY METHODS
    def appliedxoversto(self, module2):
        """
        Filters the objxovers attribute for FullXover only connecting to module2
        :param module2: Target module
        :return: list of FullXover
        """
        applied = []
        for fx in self.objxovers:
            try:
                fx.onmodule(module2)
                applied.append(fx)
            except:
                pass
        return applied

    def get_nicks(self):
        """
        Returns a list of nicks present on the module
        :return: List
        """
        nick_list = []
        for nucl in self.helix.stap:
            if nucl.toThree == -1 and nucl.__strand3__.toFive == -1:
                # nick_list.append(self.get_top_origami().get_nick(nucl))
                nuclpair = crossover.NuclPair(nucl, nucl.__strand3__, nicked=True)
                nick_list.append(nicking.Nick(nuclpair, active=False))
        return nick_list

    def get_note(self):
        """
        Returns the class Note containing module orientation and indexing
        :return: attribute self.note
        """
        return self.note


class FreeShape(Shape):
    """
    Not yet implemented.
    DNA double helix module of any shape
    Also used for compositions of piecewise modules
    """

    def __init__(self, bp, tbx, numxovers, height, dirBit):
        log.debug(__name__, "Creating FREEFORM (bp={}, "
                            "tbx={}, "
                            "xovers={}, "
                            "h={}, "
                            "dir={})".format(bp,
                                             tbx,
                                             numxovers,
                                             height,
                                             dirBit))
        Shape.__init__(self, bp, tbx, numxovers, height, dirBit)
        # additional parameters to describe a freeform shape


class LineModule(Shape):
    """
    class LineModule
    A DNA helix along a line
    """
    def __init__(self, bp, turns, dirBit, origin, normal):
        log.debug(__name__, "Creating LINE (bp={}, "
                            "turns={}, "
                            "dir={}, "
                            "origin={}, "
                            "normal={})".format(bp,
                                               turns,
                                               dirBit,
                                               origin,
                                               normal))
        Shape.__init__(self, bp, turns, dirBit, origin, normal)
        # additional parameters to describe a Line

    def __repr__(self):
        """
        Define the representation of a LineModule
        :return: string
        """
        start = "{} {} {}".format(round(self.shape.start_point[0], 2),
                                  round(self.shape.start_point[1], 2),
                                  round(self.shape.start_point[2], 2))
        end = "{} {} {}".format(round(self.shape.end_point[0], 2),
                                round(self.shape.end_point[1], 2),
                                round(self.shape.end_point[2], 2))
        dir_vec = "{} {} {}".format(round(self.shape.dir_vec[0], 2),
                                    round(self.shape.dir_vec[1], 2),
                                    round(self.shape.dir_vec[2], 2))
        return "LineModule(\n\tbp: {}, \n\ttbx: {}, \n\tnumxovers: {}, \n\tbp/t: {} \n\tstart: {}, \n\tend: {}, " \
               "\n\tdirection: {})".format(
                self.bp, self.tbx, self.numxovers, self.bpt, start, end, dir_vec)

    def __str__(self):
        start = "{} {} {}".format(round(self.shape.start_point[0], 2),
                                  round(self.shape.start_point[1], 2),
                                  round(self.shape.start_point[2], 2))
        end = "{} {} {}".format(round(self.shape.end_point[0], 2),
                                round(self.shape.end_point[1], 2),
                                round(self.shape.end_point[2], 2))
        return "LineModule [{}] [{}]".format(start, end)


class ArcRadius:
    def __get__(self, instance, owner):
        return instance.shape.radius

    def __set__(self, instance, value):
        raise AttributeError("Radius cannot be set.")


class ArcModule(Shape):
    """
    class ArcModule
    A DNA helix along an arc
    """
    # additional parameters to describe an arc
    radius = ArcRadius()

    def __init__(self, bp, turns, dirBit, origin, normal):
        log.debug(__name__, "Creating ARC (bp={}, "
                            "turns={}, "
                            "dir={}), "
                            "origin={}, "
                            "normal={}".format(bp,
                                               turns,
                                               dirBit,
                                               origin,
                                               normal))
        Shape.__init__(self, bp, turns, dirBit, origin, normal)

    def getsimpleposition(self):
        return np.array([self.helix.scaf[0].center, self.helix.scaf[-1].center])

    def __info__(self):
        """
        Show the full parameters of an ArcModule
        :return: string
        """
        center = "{} {} {}".format(round(self.shape.center[0], 2),
                                   round(self.shape.center[1], 2),
                                   round(self.shape.center[2], 2))
        return "ArcModule(\n\tbp: {}," \
               "\n\ttbx: {}," \
               "\n\tnumxovers: {}," \
               "\n\tbp/t: {}," \
               "\n\tradius: {:.2f}," \
               "\n\tcenter: {}," \
               "\n\tangle1: {}," \
               "\n\tangle2: {}".format(self.bp,
                                       self.tbx,
                                       self.numxovers,
                                       round(self.bpt, 2),
                                       round(self.shape.radius, 2),
                                       center,
                                       self.shape.init_angle,
                                       self.shape.arc_angle)

    def __repr__(self):
        """
        Show just the geometry of the module, enough distinct values to indicate which Arc it is
        :return: string
        """
        center = "{} {} {}".format(round(self.shape.center[0], 2),
                                   round(self.shape.center[1], 2),
                                   round(self.shape.center[2], 2))
        return "ArcModule [{}] [{} - {}] [{}] [{}]".format(center, self.shape.init_angle, self.shape.arc_angle,
                                                           self.bp, self.note)


class RingModule(Shape):
    """
    class RingModule
    A DNA helix formed into a circle shape

    10/08/2018
    Sending self down to DNAHelix to maintain topology
    """
    def __init__(self, bp, turns, dirBit, origin, normal, **kwargs):
        self.numxovers = kwargs.get("numxovers", None)
        self.tbx = kwargs.get("tbx", None)
        log.debug(__name__, "Creating RING (bp={}, "
                            "turns={}, "
                            "dir={}), "
                            "origin={}, "
                            "normal={}".format(bp,
                                               turns,
                                               dirBit,
                                               origin,
                                               normal))
        Shape.__init__(self, bp, turns, dirBit, origin, normal)
        # Redefine some backwards compatible values
        if self.numxovers and self.tbx:
            self.bbx = self.bp / self.numxovers  # bases between crossovers
            self.bpt = self.bbx / self.tbx  # bases per turn
            self.apb = 360.0 / self.bpt  # angle per base
        # additional properties to describe a Ring
        self.radius = mymath.bp2radius(bp)
        self.numid = -1

    def __repr__(self):
        return "RingModule [{}] [{}] [{}]".format(self.up().origin, self.bp, self.note)

    def __info__(self):
        return "RingModule(\n\tbp: {}," \
               "\n\ttbx: {}," \
               "\n\tnumxovers: {}," \
               "\n\tbp/t: {}," \
               "\n\tradius: {:.2f}," \
               "\n\tyaw: {},".format(self.bp,
                                     self.tbx,
                                     self.numxovers,
                                     round(self.bpt, 2),
                                     round(self.shape.radius, 2),
                                     self.up().yaw)

    def getposition_exact(self):
        """
        data call to get the position of the rings
        :return: array(radius, height)
        """
        return np.array((self.radius, self.height))

    def getposition_round(self):
        """
        more human readable version of getposition
        :return: array(radius, height)
        """
        return np.array((round(self.radius, 4), round(self.height, 4)))
    
    def getsimpleposition(self):
        """
        return the position as a base pair value instead of radius
        :return: array(bp, height)
        """
        return np.array((self.bp, round(self.height, 2)))

    def rotate(self, angle):
        """rotate physically moves the position of each nucleotide
        :param angle: int degrees to rotate"""
        for s in ['scaf', 'stap']:
            for nucl in self.helix.getstrand(s):        
                v = np.array([[nucl.center[0]], [nucl.center[1]]])
                v = np.transpose(rmatrix.rotate(v, angle))[0]
                duvec = np.array([[nucl.duvec[0]], [nucl.duvec[1]]])
                duvec = np.transpose(rmatrix.rotate(duvec, angle))[0]

                center = np.array([v[0], v[1], nucl.center[2]])
                duvec = np.array([duvec[0], duvec[1], nucl.duvec[2]])
                bX = mymath.unitvector(center)
                bY = self.shape.normal
                nucl.newposition(center, bX, bY, duvec)
        log.debug(__name__, "Shifted Ring({},{}) by {}°.".format(self.bp, self.height, angle))


if __name__ == "__main__":
    pass