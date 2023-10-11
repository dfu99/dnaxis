from .helper import dnaconnector, mymath, mytrig, log
from . import nucleotide, modules

import numpy as np

class Xover:
    """
    General Xover methods
    """

    def is_anti_parallel(self):
        """
        Checks if this crossover is a parallel or anti-parallel crossover
        :return: True/False
        """
        # Only need to check one half crossover
        # Dot product the direction vectors.
        # If the difference is more than 180 degrees then the strands are anti parallel
        # If the two half xovers are at a large separation distance, then check if the strands are also acyclic
        #   and if so, ignore because this is a gap-spanning crossover that will probably return a false negative
        if np.linalg.norm(self.n1.n5.center - self.n1.n3.center) > 3:  # and noncyclic.gap_exists(m1, m2):
            return True
        d1 = self.n1.n5.duvec
        d2 = self.n2.n3.duvec
        dp = np.dot(d1, d2)
        t = np.degrees(np.arccos(dp / (np.linalg.norm(d1) * np.linalg.norm(d2))))
        if t >= 90:
            return True
        else:
            # highlight_nucl(self.n1.n5)
            return False


class FullXover(Xover):
    """
    Forms a full crossover connection between pairs of nucleotide pairs
    Input two nucleotides (each is the 5' base) of two adjacent strands
    Assumes n1_5 and n2_5 are diagonal positions on adjacent anti-parallel helices
    ======53=======
    ======35=======
    ======53=======
    A full crossover is 2 half crossovers
    """

    def __init__(self, n1, n2):
        self.n1 = n1
        self.n2 = n2
        (n1_5, n1_3) = n1.to_tuple
        (n2_5, n2_3) = n2.to_tuple
        self.xovers = (HalfXover(n1_5, n2_3), HalfXover(n2_5, n1_3))

    # restores old connections
    def undo(self):
        (n1_5, n1_3) = self.n1.to_tuple
        (n2_5, n2_3) = self.n2.to_tuple
        # get the rings
        r1 = self.n1.up()
        r2 = self.n2.up()
        # remove the crossover from the rings
        r1.delxover(self)
        r2.delxover(self)
        # remove reference to HalfXovers
        self.xovers = None
        # reconnect original connections
        log.debug(__name__, "Undo crossover", self.n1.numid, "-", self.n2.numid)
        dnaconnector.nuclconnect(n1_5, n1_3)
        dnaconnector.nuclconnect(n2_5, n2_3)

    def onmodule(self, module):
        """
        The crossover does not have any built-in directionality, thus we do not know which pair of nucleotides
            is on which side of the crossover (and on which strand)
        Helps discern which NuclPair is on which module
        :param module: the module object
        :return: the NuclPair that is on module
        """
        m1 = self.n1.up()
        m2 = self.n2.up()
        if m1 == module:
            return self.n1
        elif m2 == module:
            return self.n2
        else:
            raise ValueError("This crossover does not exist on Ring({},{}).".format(module.bp, module.height))

    def __contains__(self, value):
        """
        Checks if Nucleotide value is present in the XoverPair
        :param value:
        :return:
        """
        if value in (self.n1.n5, self.n1.n3, self.n2.n5, self.n2.n3):
            return True
        else:
            return False

    def spans(self, m1, m2):
        """
        Returns True/False whether crossover spans the input strands
        :param m1: Any shape or strand object
        :param m2: Any shape or strand object
        :return: True/False
        """
        if issubclass(type(m1), modules.Shape) and issubclass(type(m2), modules.Shape):
            nt1 = self.n1.n5.get_top_module()
            nt2 = self.n2.n5.get_top_module()

        elif type(m1) == modules.JoinedStrand and type(m2) == modules.JoinedStrand:
            nt1 = self.n1.n5.get_top_strand()
            nt2 = self.n2.n5.get_top_strand()

        # Checks equality by their topology value (should be correctly set)
        if all(x in [m1, m2] for x in [nt1, nt2]):
            return True
        else:
            return False

    def has_nucl(self, nucl):
        """
        Returns True/False whether Nucleotide nucl comprises the FullXover or not
        :param nucl: Nucleotide
        :return: True/False
        """
        if nucl in [self.n1.n5, self.n1.n3, self.n2.n5, self.n2.n3]:
            return True
        else:
            return False

    def settype(self, choose_strand):
        if choose_strand not in ['scaf', 'stap']:
            raise ValueError("Strand type must be 'scaf' or 'stap'.")
        self.strandtype = choose_strand

    def __repr__(self):
        h1 = self.xovers[0]
        h2 = self.xovers[1]
        return "\n{:^15} {:^15}\n{:^15} {:^15}\n{:^15} {:^15}".format(
            h1.nucls[0].numid, h2.nucls[1].numid, "|", "|", h1.nucls[1].numid, h2.nucls[0].numid)

    def __eq__(self, other):
        """
        Checks equivalence by matching all 4 nucleotide crossover points
        :param other: FullXover
        :return: True/False
        """
        if type(other) != type(self):
            return False
        x1n1n5 = self.n1.n5
        x1n1n3 = self.n1.n3
        x1n2n5 = self.n2.n5
        x1n2n3 = self.n2.n3
        x1 = [x1n1n5, x1n1n3, x1n2n5, x1n2n3]
        x2n1n5 = other.n1.n5
        x2n1n3 = other.n1.n3
        x2n2n5 = other.n2.n5
        x2n2n3 = other.n2.n3
        x2 = [x2n1n5, x2n1n3, x2n2n5, x2n2n3]
        if all(nid in x2 for nid in x1) and all(nid in x1 for nid in x2):
            return True
        else:
            return False

    def complement(self):
        """
        Returns a virtual complement as an XoverPair
        :return: XoverPair
        """
        (n1_5, n1_3) = self.n1.to_tuple
        (n2_5, n2_3) = self.n2.to_tuple
        n1 = NuclPair(n1_3.Comp, n1_5.Comp, nicked=True)
        n2 = NuclPair(n2_3.Comp, n2_5.Comp, nicked=True)
        xp = XoverPair(n1, n2)
        return xp


class HalfXover:
    def __init__(self, n_5, n_3):
        self.nucls = (n_5, n_3)
        dnaconnector.nuclconnect(n_5, n_3)

    def __repr__(self):
        return "{}-{}".format(self.nucls[0].numid, self.nucls[1].numid)


class XoverPair(Xover):
    """
    facilitate the format of an xover pair as type tuple
    define: an xover pair is a pair of nucleotide pairs on adjacent helices
      where an xover could be placed
    an xover pair should be formatted:
    (<nucleotide pair (5',3') #1>,
    <nucleotide pair (5',3') #2>,
    <distance between the nucleotide pairs>,
    <angular position of nucleotide pair #1>,
    <angular position of nucleotide pair #2>)
    """

    def __init__(self, n1, n2):
        self.n1 = n1  # n1 is NuclPair object
        self.n2 = n2  # n2 is NuclPair object
        self.dist = abs(mytrig.angle_diff(self.n1.pos, self.n2.pos))

    def fget(self):
        return self.n1, self.n2, self.dist, self.n1.pos, self.n2.pos

    def fset(self, *value):
        (self.n1, self.n2, self.dist) = value

    def apply(self):
        return FullXover(self.n1, self.n2)

    def get_halfxover(self, nucl):
        """
        Gets the HalfXover nucleotides corresponding to the input Nucleotide
        :param nucl: Nucleotide object
        :return: tuple, in order of found nucl first
        """
        (n1_5, n1_3) = self.n1.to_tuple
        (n2_5, n2_3) = self.n2.to_tuple

        if nucl == n1_5:
            return n1_5, n2_3
        elif nucl == n2_3:
            return n2_3, n1_5
        elif nucl == n2_5:
            return n2_5, n1_3
        else:
            return n1_3, n2_5

    def onmodule(self, module):
        """
        The crossover does not have any built-in directionality, thus we do not know which pair of nucleotides
            is on which side of the crossover (and on which strand)
        Helps discern which NuclPair is on which module
        :param module: the module object
        :return: the NuclPair that is on module
        """
        m1 = self.n1.up()
        m2 = self.n2.up()
        if m1 == module:
            return self.n1
        elif m2 == module:
            return self.n2
        else:
            raise ValueError("This crossover does not exist on Ring({},{}).".format(module.bp, module.height))

    def __eq__(self, other):
        """
        Checks equivalence by matching all 4 nucleotide crossover points
        :param other: XoverPair
        :return: True/False
        """
        if type(other) != type(self):
            return False
        x1n1n5 = self.n1.n5
        x1n1n3 = self.n1.n3
        x1n2n5 = self.n2.n5
        x1n2n3 = self.n2.n3
        x1 = [x1n1n5, x1n1n3, x1n2n5, x1n2n3]
        x2n1n5 = other.n1.n5
        x2n1n3 = other.n1.n3
        x2n2n5 = other.n2.n5
        x2n2n3 = other.n2.n3
        x2 = [x2n1n5, x2n1n3, x2n2n5, x2n2n3]
        if all(nid in x2 for nid in x1) and all(nid in x1 for nid in x2):
            return True
        else:
            return False

    def get_nuclpair(self, nucl):
        """
        Gets the NuclPair nucleotides corresponding to the input Nucleotide
        :param nucl: Nucleotide object
        :return: tuple
        """
        (n1_5, n1_3) = self.n1.to_tuple
        (n2_5, n2_3) = self.n2.to_tuple

        if nucl == n1_5:
            return n1_5, n1_3
        elif nucl == n1_3:
            return n1_3, n1_5
        elif nucl == n2_3:
            return n2_3, n2_5
        else:
            return n2_5, n2_3

    def get_nucls(self):
        return self.n1.n5, self.n1.n3, self.n2.n5, self.n2.n3

    def __contains__(self, value):
        """
        Checks if Nucleotide value is present in the XoverPair
        :param value:
        :return:
        """
        if value in (self.n1.n5, self.n1.n3, self.n2.n5, self.n2.n3):
            return True
        else:
            return False

    def spans(self, m1, m2):
        """
        Returns True/False whether crossover spans the input modules
        :param m1: modules.Shape
        :param m2: modules.Shape
        :return: True/False
        """
        if issubclass(type(m1), modules.Shape) and issubclass(type(m2), modules.Shape):
            nt1 = self.n1.n5.get_top_module()
            nt2 = self.n2.n5.get_top_module()

        elif type(m1) == modules.JoinedStrand and type(m2) == modules.JoinedStrand:
            nt1 = self.n1.n5.get_top_strand()
            nt2 = self.n2.n5.get_top_strand()

        # Checks equality by their topology value (should be correctly set)
        if all(x in [m1, m2] for x in [nt1, nt2]):
            return True
        else:
            return False

    def __repr__(self):
        return "\n{:^15} {:^15} {}\n{:^15} {:^15}\n{:^15} {:^15} {} Bond Length: {}".format(
            self.n1.n5.numid, self.n1.n3.numid, self.n1.n3.get_top_module(), "|", "|", self.n2.n3.numid,
            self.n2.n5.numid, self.n2.n5.get_top_module(), round(self.bondlen(), 2))

    def bondlen(self):
        return np.linalg.norm(self.n2.coords - self.n1.coords)

    to_tuple = property(fget, fset)


class NuclPairPos:
    """
    Returns the angular position of a NuclPair
    ONLY FOR RING MODULES
    """

    def __get__(self, instance, owner):
        return mymath.getmdpt(mymath.nucl2angle(instance.n5), mymath.nucl2angle(instance.n3))

    def __set__(self, instance, value):
        raise AttributeError("Position cannot be set.")


class NuclPairTheta:
    """
    Returns the average theta angle value of two adjacent nucleotides of a NuclPair
    """

    def __get__(self, instance, owner):
        return mymath.getmdpt(instance.n5.theta, instance.n3.theta)

    def __set__(self, instance, value):
        raise AttributeError("Theta cannot be set.")


class NuclPairCoords:
    """
    Returns the midpoint of the centers of nucleotides in a NuclPair
    """

    def __get__(self, instance, owner):
        return (instance.n5.center + instance.n3.center) / 2

    def __set__(self, instance, value):
        raise AttributeError("Coordinates cannot be set.")


class NuclPairBisectingVector:
    """
    Returns the bisecting vector of two nucleotides representing their average direction
    """

    def __get__(self, instance, owner):
        a = instance.n5.nucl_vec
        b = instance.n3.nucl_vec
        c = np.linalg.norm(a) * b + np.linalg.norm(b) * a
        return c

    def __set__(self, instance, value):
        raise AttributeError("BisectingVector cannot be set.")


class NuclPair:
    """
    Facilitate the format of a nucleotide pair as type tuple
    Define: a nucleotide pair is a pair of nucleotides that are adjacent
      on the same strand
    A nucleotide pair should be formatted:
    (<5' nucleotide>, <3' nucleotide>)

    2/17/2020: Removed condition that nucleotides have to be on same Module
                Modified to being on same Plane
    """
    pos = NuclPairPos()
    theta = NuclPairTheta()
    coords = NuclPairCoords()
    vector = NuclPairBisectingVector()

    def __init__(self, *nucl, nicked=False):
        self.nicked = nicked
        if self.chkformat(*nucl):
            (self.n5, self.n3) = nucl
            self.top = self.n5.up().up()  # to ringmodule

    def fget(self):
        """
        Gets NuclTuple
        :return: 5' and 3' Nucleotide objects
        """
        return self.n5, self.n3

    def fset(self, *value):
        """
        Sets NuclTuple
        :param value: 5' and 3' Nucleotide arguments
        :return: None
        """
        if self.chkformat(self, *value):
            (self.n5, self.n3) = value

    def nget(self):
        return self.n5.numid, self.n3.numid

    def nset(self):
        pass

    # returns ring
    def up(self):
        return self.top

    def chkformat(self, *value):
        """
        Checks if the input are Nucleotide objects
        :param value: Nucleotide objects
        :param nicked: True/False for adjacency/connectivity
        :return:
        """
        # Check that input is a 2-tuple
        try:
            (nucl5, nucl3) = value
        except:
            raise ValueError("{} does not match (5',3') nucleotide format.".format(value))

        # Check that both elements are Nucleotide objects
        if not (isinstance(nucl5, nucleotide.Nucleotide) and isinstance(nucl3, nucleotide.Nucleotide)):
            raise TypeError("Cannot set {} or {} to type 'Nucleotide'".format(type(nucl3), type(nucl5)))

        # Check for adjacency or connectivity
        if self.nicked:
            if nucl5.__strand3__ != nucl3 and nucl3.__strand5__ != nucl5:
                raise RuntimeError("Nucleotides must be adjacent to form a NuclPair.")
        else:
            if nucl5.toThree != nucl3 and nucl3.toFive != nucl5:
                raise RuntimeError("Nucleotides must be connected to form a NuclPair.")

        # Check that they are on the same strand
        if nucl5.get_top_strand() != nucl3.get_top_strand():
            raise AttributeError("Nucleotides do not belong to same strand. n5 [{}] n3 [{}]".format(
                nucl5.up().up(),
                nucl3.up().up()))

        return True

    to_tuple = property(fget, fset)
    numid = property(nget, nset)


class GapNuclPair(NuclPair):
    """
    Inherits from NuclPair for crossovers spanning gaps (half crossover placed on each side)
    """

    def __init__(self, *nucl):
        if self.chkformat(*nucl):
            (self.n5, self.n3) = nucl

        # # Break off the excess helix
        # # Scaffold
        # dnaconnector.nuclbreak(self.n5)
        # dnaconnector.nuclbreak(self.n3.__strand5__)
        # # Staple
        # dnaconnector.nuclbreak(self.n5.Comp.__strand5__)
        # dnaconnector.nuclbreak(self.n3.Comp)
        # # The strand's original absolute connections also need to be terminated
        # # Scaffold
        # dnaconnector.truebreak(self.n5)
        # dnaconnector.truebreak(self.n3.__strand5__)
        # # Staple
        # dnaconnector.truebreak(self.n5.Comp.__strand5__)
        # dnaconnector.truebreak(self.n3.Comp)

        # Emulate a cycle by spanning the gap for the staple and scaffold
        dnaconnector.trueconnect(self.n5, self.n3)
        dnaconnector.trueconnect(self.n3.Comp, self.n5.Comp)
        # Reconnect scaffolds (they will be disconnected later anyways to make the crossover)
        dnaconnector.nuclconnect(self.n5, self.n3)
        # Do not reconnect staples
        # TODO: Get rid of this, this is horrible
        dnaconnector.nuclconnect(self.n3.Comp, self.n5.Comp)

        # Removed constraint for gaps
        # if self.n5.toThree != self.n3 and self.n3.toFive != self.n5:
        #     raise RuntimeError("Nucleotides must be connected to form a NuclPair.")

        if self.n5.get_top_strand() != self.n3.get_top_strand():
            raise AttributeError("Nucleotides do not belong to same strand. n5 [{}] n3 [{}]".format(
                self.n5.up().up(),
                self.n3.up().up()))
        else:
            self.top = self.n5.up().up()  # to ringmodule

    def chkformat(self, *value):
        """
        Checks if the input are Nucleotide objects
        :param value:
        :return:
        """
        try:
            (nucl5, nucl3) = value
        except:
            raise ValueError("{} does not match (5',3') nucleotide format.".format(value))
        if not (isinstance(nucl5, nucleotide.Nucleotide) and isinstance(nucl3, nucleotide.Nucleotide)):
            raise TypeError("Cannot set {} or {} to type 'Nucleotide'".format(type(nucl3), type(nucl5)))
        # Removed constraint for gaps
        # if nucl5.toThree != nucl3 or nucl3.toFive != nucl5:
        #     raise ValueError("Nucleotides {} and {} are not connected or adjacent.".format(nucl5.numid, nucl3.numid))
        # else:
        #     return True
        return True