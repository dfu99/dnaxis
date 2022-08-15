"""
base object

Author: Dan Fu
Last changed: 5/16/2017

CHANGELOG:
    4/1/2017: File created
        Initialize function created
    5/16/2017: Bug fix
        Variable ctheta should be written once with angle() function too
        Changed:
            self.ctheta = (ctheta + twist)%360.0000
            to
            self.ctheta = (ctheta + twist)%360.0000
    8/14/2017
        Removed ctheta
    8/23/2017
        Added nucl variable tracking base letter (ATCG)
    3/30/2018
        Move calculation of theta and twist from angle() function to scaffoldCW/CCW() functions in dnacircle.py. This should be assignment only.
    4/2/2018:
        Rename to nucleotide.py and class to InitNucl
        Add x,y,z,theta as initialization parameters
        Add oxDNA values in initialization parameters
        Change toThree, toFive, Comp to initialize as -1 instead of 0
        Remove r (radius)
        Remove hexid because it is the same as numid
        Remove angle function (no longer necessary if theta is initialized)
    8/2/2018:
        Bugfix: Theta should always be mod 360
    8/19/2018:
        Add a printid parameter to faciliate printout to oxDNA format
    9/12/2018:
        Renamed class to Nucleotide
    10/08/2018
        Added topology argument and instance variable
    10/13/2018
        Broken into smaller methods, setoxdna, setxyz, settop, initrandnucl
        self.update added to init
        theta, direction added to init args
        [nucl_vec, nucl_pos] renamed to [_nucl_vec, nucl_vec]
"""
import random
import math
import numpy as np
from .helper import mymath, log
from app import config


# returns the complement of a nucleotide
def comp_nucl(n):
    if n.upper() == 'A':
        return 'T'
    elif n.upper() == 'T':
        return 'A'
    elif n.upper() == 'G':
        return 'C'
    elif n.upper() == 'C':
        return 'G'
    else:
        log.out(__name__, "ERROR: Bad input to comp_nucl.")
        raise RuntimeError("ERROR: Bad input to comp_nucl.")
        return 'X'


def randnucl():
    """
    Gets random nucleotide
    :return: A T C G
    """
    nucl_bank = ['A', 'T', 'C', 'G']
    random.shuffle(nucl_bank)
    return nucl_bank[0]


class AsymTheta:
    # DEPRECATED
    def __get__(self, instance, owner):
        # (x, y, z) = (instance.x, instance.y, instance.z)
        # (x, y, z) = np.array([x, y, z])-instance.center
        (x, y, z) = instance.nucl_vec
        (rho, theta, phi) = mymath.cart2sph(x, y, z)
        phi = math.degrees(phi)
        return phi

    def __set__(self, instance, value):
        raise AttributeError("Asymmetric theta cannot be set.")


class Nucleotide(object):
    asymtheta = AsymTheta()  # DEPRECATED

    def __init__(self, numid, node, bX, bY, duvec, theta, direction, top):
        # topology
        self.settop(top)
        # base
        self.initrandnucl()
        # ID
        self.numid = numid
        # connection info
        self.toThree = -1
        self.toFive = -1
        self.Comp = -1

        self.node = node  # node preserves the shape connection graph
        center = node.coordinates
        (self.rho, phi) = mymath.cart2pol(center[0], center[1])
        self.phi = math.degrees(phi) % 360
        
        # strand connection info
        self.__strand3__ = -1
        self.__strand5__ = -1

        # oxDNA parameters
        self.setoxdna(center)
        self.printid = -1  # helps for printing to top file

        # basis vectors for nucleotides along strand
        self.bX = bX
        self.bY = bY
        self.duvec = duvec  # direction unit vector

        self.update(theta, direction)
        
        # caDNAno parameters
        self.cdna_num = -1
        self.cdna_ind = -1

        # Piecewise strand topology
        self.strand = None

    def get_top_strand(self):
        """
        :return: Containing JoinedStrand
        """
        return self.strand

    def get_top_helix(self):
        """
        :return: Containing DNAHelix
        """
        return self.up()

    def get_top_module(self):
        """
        :return: Containing ShapeModule
        """
        return self.get_top_helix().up()

    def get_top_plane(self):
        """
        :return: Containing Plane
        """
        return self.get_top_module().up()

    def get_top_origami(self):
        """
        :return: Containing Origami
        """
        return self.get_top_plane().up()

    def update(self, theta, direction):
        """calculates the position of the nucleotide in x, y, z and oxdna formats for a new angle position theta
        :param theta: int [0,360] degrees new angle of nucleotide w.r.t. same center
        :param direction: only used for setting direction during __init__
        :return: None"""
        bdna_r = config.BDNA_R
        # Where on the circular slice is the base? 
        # 0 degree corresponds to 3 o clock.
        # + CCW, - CW
        self.theta = theta % 360
        self.direction = direction
        # x and y magnitudes on flat plane
        _nucl_vec = [bdna_r*math.cos(math.radians(theta)), bdna_r*math.sin(math.radians(theta))]
        # project on direction
        nucl_vec = _nucl_vec[0]*self.bX+_nucl_vec[1]*self.bY
        self.nucl_vec = nucl_vec
        # calculate spatial data
        [x, y, z] = nucl_vec+self.center
        self.setxyz(x, y, z)
        
        # calculate oxDNA data
        p = self.center+nucl_vec*0.6
        bbv = -nucl_vec
        nv = -self.duvec
        self.setoxdna(self.center, p, bbv, nv)
        
    def newposition(self, center, bX, bY, duvec):
        """calculates the position of the nucleotide in x, y, z and oxdna formats for a new center position
        :param center: x, y, z array
        :param bX: new X basis vector
        :param bY: new Y basis vector
        :param duvec: new 3D direction unit vector
        :return: None"""
        self.setoxdna(center)
        self.bX = bX
        self.bY = bY
        self.duvec = duvec
        self.update(self.theta, self.direction)
    
    def setoxdna(self, c=None, p=None, bbv=None, nv=None):
        self.center = c
        self.position = p
        self.backbonevec = bbv
        self.normvec = nv

    def setxyz(self, x=None, y=None, z=None):
        self.x = x
        self.y = y
        self.z = z
        
    def settop(self, top):
        self.top = top
        
    def initrandnucl(self):
        """
        Initializes nucl base label
        """
        self.nucl = randnucl()
        
    def up(self):
        """
        up one level to higher class
        DNAHelix object contains Nucleotide
        :return: DNAHelix address
        """
        return self.top

    def __repr__(self):
        return str(self.numid)

    def is_on_gap(self):
        """
        True if Nucleotide is adjacent to a gap, False otherwise
        :return:
        """
        if np.linalg.norm(self.center - self.__strand5__.center) > 3 or \
                np.linalg.norm(self.center - self.__strand3__.center) > 3:
            return True
        return False

    def strand_type(self):
        """
        :return: Returns 'scaf' or 'stap'
        """
        if self in self.get_top_helix().scaf:
            return 'scaf'
        elif self in self.get_top_helix().stap:
            return 'stap'
        else:
            log.system("E0x2")
            log.system(self.numid)
            log.system(self.get_top_module())
            log.system("RuntimeError: This nucleotide is not linked properly to a DNAHelix")
            raise RuntimeError("This nucleotide is not linked properly to a DNAHelix.")
