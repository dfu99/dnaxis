# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 19:52:31 2018

@author: Dan
DNAHelix.py populates points with nucleotides

4/2/2018:
    Split previous dnacircle.py as two objects:
        One object initiates the shape
        Another object populates the bases
    Primary difference is to split shape formation from nucleotides (i.e. makecircle, scaffoldCW/CCW, etc.)
    This is expected to be more general and more concise.
    Can take any path as input, i.e circles, lines, or even freeform.
    Instead of IDs indicating 5', 3', and base paired bases, use the nucleotide object
    
7/16/2018:
    - Bugfix:
        input to InitNucl is (curr_base_index,bX,bY,udir_vec,this_node)
        but __init__ is actually (self,numid,center,bX,bY,duvec)
9/12/2018:
    - Renamed InitNucl to Nucleotide
10/08/2018
    Split up the DNAHelix constructor into additional methods
    Added topology argument and instance variable
"""
import math

import numpy as np

from . import nucleotide, makeshape
from app import config
from .helper import mymath, mytrig

class Strand:
    def __init__(self, nucl_list):
        self.bases = nucl_list


class DNAHelix(object):
    # constructor defines initialization values
    # and creates strand lists for scaffold and staple
    def __init__(self, shape_object, angle_per_base, init_angle, init_base_index, direction, top):
        """
        Initializes the DNA helix according to the input shape and helix geometry parameters
        :param shape_object: shape node graph
        :param angle_per_base: angular offset between adjacent nucleotides
        :param init_angle: angle of first nucleotide (possibly redundant)
        :param init_base_index: nucleotide index numbering to start from
        :param direction: 5'-3' direction indicator. Always uses SCAFFOLD strand as reference
        :param top: module containing the helix
        """
        # Strand variables
        # @ Invariant: Strand sequence list always begins with 5' nucleotide
        self.scaf = []
        self.stap = []
        
        # instance variable constants for object
        self.angle_per_base = angle_per_base
        self.init_angle = init_angle
        self.graph = shape_object.graph
        self.direction = direction
        self.num_points = len(self.graph)
        self.top = top
        self.normal = shape_object.normal
        
        # Initialize counters
        sc_theta = init_angle
        curr_base_index = init_base_index

        # 1. Derive sub graphs
        if direction:  # True, follow the graph nodes
            sc_nodegraph = self.graph
            st_nodegraph = self.graph[::-1]  # Always has complementary 5'-3' direction to scaffold
        else:  # False, reverse the graph nodes
            sc_nodegraph = self.graph[::-1]
            st_nodegraph = self.graph

        # 2. Populate nucleotides
        # First, make the scaffold strand
        # At each center point stated by coordinates, treat possible positions for nucleotides in a circle
        # existing on a plane normal to its direction vector (dir_vec)
        for node in sc_nodegraph:  # do not iterate past last coordinate
            # Vector from this point to next point
            try:
                if direction:
                    next_node = node.next_point
                else:
                    next_node = node.prev_point
                next_coords = next_node.coordinates
            except AttributeError:  # If nodes are not cyclic, last node will be an 'int'
                # Assume coordinates continue in the same direction to make the last unit direction vector (udir_vec)
                next_coords = this_coords + 2 * udir_vec * config.AXIAL_RISE
            this_node = node
            this_coords = this_node.coordinates
            if np.array_equal(this_coords, next_coords):  # Error handler
                raise RuntimeError("Error: The node did not proceed. "
                                   "Next coordinates {}, "
                                   "current coordinates: {}, "
                                   "direction vector: {}".format(
                                    next_coords, this_coords, udir_vec))
            # Nucleotide plane basis vectors
            udir_vec = mymath.unitvector(next_coords-this_coords)
            (bX, bY) = mymath.get_basis_vec(shape_object.normal, udir_vec)
            # Initialize scaffold nucleotides
            self.populate(self.scaf,
                          curr_base_index,
                          sc_theta,
                          this_node,
                          bX, bY, udir_vec, self.direction)
            # Increment the index and angle
            curr_base_index += 1
            sc_theta = (sc_theta-angle_per_base) % 360.0

        # After scaffold nucleotides are populated, do the same for staple nucleotides
        for node in st_nodegraph:  # Do not iterate past last coordinate
            # Vector from this point to next point
            try:
                if direction:
                    next_node = node.prev_point
                else:
                    next_node = node.next_point
                next_coords = next_node.coordinates
            except AttributeError:  # If nodes are not cyclic, last node will be an 'int'
                # Assume coordinates continue in the same direction to make the last unit direction vector (udir_vec)
                next_coords = this_coords + 2 * udir_vec * config.AXIAL_RISE
            this_node = node
            this_coords = this_node.coordinates
            if np.array_equal(this_coords, next_coords):  # Error handler
                raise RuntimeError("Error: The node did not proceed. "
                                   "Next coordinates {}, "
                                   "current coordinates: {}, "
                                   "direction vector: {}".format(
                                    next_coords, this_coords, udir_vec))
            # Nucleotide basis vectors
            udir_vec = mymath.unitvector(next_coords - this_coords)
            (bX, bY) = mymath.get_basis_vec(shape_object.normal, udir_vec)
            # Get the staple theta relevant to the scaffold nucleotide.
            # One strand always lags the other
            st_theta = mymath.lh2rh2((this_node.scaf.theta + config.NUCLLAGANGLE) % 360)
            # Initialize staple nucleotides
            self.populate(self.stap,
                          curr_base_index,
                          st_theta,
                          this_node,
                          bX, bY, udir_vec, not self.direction)
            # Increment the index
            curr_base_index += 1

        # 3. Traverse node graph and set Complementary Nucleotides
        for node in self.graph:
            node.scaf.Comp = node.stap
            node.stap.Comp = node.scaf

        # 4. Traverse the sub graphs and set 5' and 3' adjacency
        # For the scaffold graph
        for node in sc_nodegraph:
            try:
                if direction:
                    node.scaf.toThree = node.next_point.scaf
                    node.scaf.__strand3__ = node.next_point.scaf
                else:
                    node.scaf.toThree = node.prev_point.scaf
                    node.scaf.__strand3__ = node.prev_point.scaf
            except AttributeError:
                node.scaf.toThree = -1
                node.scaf.__strand3__ = -1
            try:
                if direction:
                    node.scaf.toFive = node.prev_point.scaf
                    node.scaf.__strand5__ = node.prev_point.scaf
                else:
                    node.scaf.toFive = node.next_point.scaf
                    node.scaf.__strand5__ = node.next_point.scaf
            except AttributeError:
                node.scaf.toFive = -1
                node.scaf.__strand5__ = -1

        # For the staple graph
        for node in st_nodegraph:
            try:
                if direction:
                    node.stap.toThree = node.prev_point.stap
                    node.stap.__strand3__ = node.prev_point.stap
                else:
                    node.stap.toThree = node.next_point.stap
                    node.stap.__strand3__ = node.next_point.stap
            except AttributeError:
                node.stap.toThree = -1
                node.stap.__strand3__ = -1
            try:
                if direction:
                    node.stap.toFive = node.next_point.stap
                    node.stap.__strand5__ = node.next_point.stap
                else:
                    node.stap.toFive = node.prev_point.stap
                    node.stap.__strand5__ = node.prev_point.stap
            except AttributeError:
                node.stap.toFive = -1
                node.stap.__strand5__ = -1

    def populate(self, strand, ind, theta, this_node, bX, bY, udir_vec, direction):
        """
        Add nucleotide objects to strands
        Update positions of nucleotides
        :param strand: This DNAHelix's scaf or stap attribute
        :param ind: int Current index number
        :param theta: float Angle
        :param this_node: Point1d Graph node
        :param bX: 3-array vector X Basis for nucleotide positioning
        :param bY: 3-array vector Y Basis for nucleotide positioning
        :param udir_vec: 3-array vector Normal direction to next nucleotide
        :param direction: binary Corresponds to 5' or 3' direction of the strand
        :return: None
        """
        new_nucl = nucleotide.Nucleotide(ind, this_node, bX, bY, udir_vec, theta, direction, self)
        strand.append(new_nucl)
        if strand is self.stap:
            this_node.stap = new_nucl  # Link to Point1d
            if this_node.stap == 'ERROR':  # Error flagging if node pairing to nucleotide was not set
                raise ValueError("Connection error!")
        else:
            this_node.scaf = new_nucl  # Link to Point1d
            if this_node.scaf == 'ERROR':  # Error flagging if node pairing to nucleotide was not set
                raise ValueError("Connection error!")

    def up(self):
        '''
        up one level to higher level containing object
        Ring object contains DNAHelix
        :return: Ring address
        '''
        return self.top

    # shifts each nucleotide of the entire ring by a certain angle
    def twist(self, ref_angle, ref_nucl):
        """
        Twists the entire strand to satisfy the angle position of the reference nucleotide
        :param ref_angle: int degrees new angle of reference nucleotide
        :param ref_nucl: reference nucleotide on strand
        """

        # validate that the nucl is on the helix
        if ref_nucl not in self.scaf and ref_nucl not in self.stap:
            raise ValueError("Reference nucleotide was not found.")

        nucl_angle = ref_nucl.theta
        angle_shift = mytrig.angle_diff(ref_angle, nucl_angle)

        # (Deprecated) Strand carryover method which rotates already connected strands,
        #   even those not within the current module

        # strand = strandnav.getstrand(ref_nucl)
        # sc_theta = (strand[0].theta + angle_shift) % 360
        # for scnucl in strand:
        #     scnucl.update(sc_theta, self.direction)
        #     st_theta = sc_theta-(math.pow(-1, self.direction)*config.NUCLLAGANGLE)
        #     stnucl = scnucl.Comp
        #     stnucl.update(st_theta, not self.direction)
        #     sc_theta = (sc_theta - self.angle_per_base) % 360.0
        # Module-only method
        sc_theta = (self.scaf[0].theta + angle_shift) % 360
        for scnucl in self.scaf:
            scnucl.update(sc_theta, self.direction)
            st_theta = mymath.lh2rh2((scnucl.theta + config.NUCLLAGANGLE) % 360)
            stnucl = scnucl.node.stap
            stnucl.update(st_theta, not self.direction)
            sc_theta = (sc_theta-self.angle_per_base) % 360.0

        # check
        rnd_nucl_theta = round(ref_nucl.theta, 2)
        rnd_ref_theta = round(ref_angle, 2)
        if rnd_nucl_theta != rnd_ref_theta:
            raise RuntimeError("Twist failed. The reference nucleotide {} angle {} "
                               "did not match the reference angle {} after twisting.".format(
                                ref_nucl.numid, ref_nucl.theta, ref_angle))
    
    def getstrand(self, choose_strand):
        if choose_strand == 'scaf':
            return self.scaf
        elif choose_strand == 'stap':
            return self.stap
        else:
            raise ValueError("Strand choice '{}' does not exist. Expected 'scaf' or 'stap' strand.")

    def insertNucl(self, nucl1, nucl2):
        """
        insert adds a nucleotide on the helix (both scaf and stap) between nucleotides given by nucl1 and nucl2
        :param nucl1: type Nucleotide
        :param nucl2: type Nucleotide
        :return: None
        """
        # check if nucl1 and nucl2 are continuous
        if not (nucl1.toThree or nucl1.toFive == nucl2):
            raise ValueError("Nucleotides are not continuous.")

        # check which strand it's on
        if nucl1 in self.stap:
            nucl1 = nucl1.Comp
            nucl2 = nucl2.Comp
        # find the lesser index within the scaffold list
        ind1 = self.scaf.index(nucl1)
        ind2 = self.scaf.index(nucl2)
        if ind1<ind2:
            ind = ind1
        else:
            ind = ind2
        # determine which nucl are the 5' and 3' side
        if nucl1.toThree == nucl2:
            nucl5 = nucl1
            nucl3 = nucl2
        elif nucl2.toThree == nucl1:
            nucl5 = nucl2
            nucl3 = nucl1
        else:
            raise ValueError("These nucleotides are not connected!")
        # set the numid
        ins_id = nucl5.numid+100000
        # find the new center
        # convert nucl5 and nucl3 into polar to find their angle
        (rho1, phi1) = mymath.cart2pol(nucl5.center[0], nucl5.center[1])
        (rho2, phi2) = mymath.cart2pol(nucl3.center[0], nucl3.center[1])
        phi1 = math.degrees(phi1) % 360
        phi2 = math.degrees(phi2) % 360
        # find the bisecting angle
        ins_phi = math.radians(mymath.getmdpt(phi1, phi2))
        # keep the same radius. average in case of some variance
        ins_rho = (rho1+rho2)/2
        # get the new center
        (ins_x, ins_y) = mymath.pol2cart(ins_rho, ins_phi)
        ins_z = nucl5.center[2]

        ins_center = np.array([ins_x, ins_y, ins_z])
        ins_node = makeshape.Point1d(-1, -1, ins_center)
        # get the new direction unit vector to the 3' nucl
        ins_duvec = mymath.unitvector(nucl3.center - ins_center)
        (ins_bX, ins_bY) = mymath.get_basis_vec(self.normal, ins_duvec)
        # calculate the theta position of the nucleobase for the scaffold strand
        ins_sctheta = mymath.getmdpt(nucl5.theta, nucl3.theta)
        # create the nucleotide and insert it into the strand
        ins_scnucl = nucleotide.Nucleotide(ins_id,
                                           ins_node, ins_bX, ins_bY, ins_duvec, ins_sctheta, self.direction, self)
        self.scaf.insert(ind, ins_scnucl)
        # calculate the theta position of the staple nucleobase w.r.t. the scaffold nucleobase
        ins_sttheta = ins_sctheta-(math.pow(-1, self.direction)*config.NUCLLAGANGLE)
        # create the nucleotide and insert it into the strand
        ins_stnucl = nucleotide.Nucleotide(ins_id+100000,
                                           ins_node, ins_bX, ins_bY, ins_duvec, ins_sttheta, not self.direction, self)
        self.stap.insert(ind, ins_stnucl)
        # fix the connections
        ins_scnucl.Comp = ins_stnucl
        ins_stnucl.Comp = ins_scnucl

        ins_scnucl.toFive = nucl5
        nucl5.toThree = ins_scnucl
        ins_scnucl.toThree = nucl3
        nucl3.toFive = ins_scnucl

        comp3 = nucl5.Comp
        comp5 = nucl3.Comp

        ins_stnucl.toFive = comp5
        comp5.toThree = ins_stnucl
        ins_stnucl.toThree = comp3
        comp3.toFive = ins_stnucl

    def deleteNucl(self, nucl):
        """
        delete removes a nucleotide on the helix (both scaf and stap) at nucl
        :param strand: scaf or stap list attribute of DNAHelix
        :param nucl: type Nucleotide
        :return: None
        """
        # check if nucl is on scaffold
        if nucl in self.scaf:
            nucl5 = nucl.toFive
            nucl3 = nucl.toThree
            self.scaf.remove(nucl)
            nucl5.toThree = nucl3
            nucl3.toFive = nucl5

            comp = nucl.Comp
            comp5 = comp.toFive
            comp3 = comp.toThree
            self.stap.remove(comp)
            try:
                comp5.toThree = comp3
            except:
                pass
            try:
                comp3.toFive = comp5
            except:
                pass
        # check if nucl is on staple
        elif nucl in self.stap:
            nucl5 = nucl.toFive
            nucl3 = nucl.toThree
            self.stap.remove(nucl)
            nucl5.toThree = nucl3
            nucl3.toFive = nucl5

            comp = nucl.Comp
            comp5 = comp.toFive
            comp3 = comp.toThree
            self.scaf.remove(comp)
            try:
                comp5.toThree = comp3
            except:
                pass
            try:
                comp3.toFive = comp5
            except:
                pass
        else:
            raise ValueError("Nucleotide is not found on helix.")

