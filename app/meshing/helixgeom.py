# scaleshape.m expands a shape via polar coordinates according to the
# scalingFactor

import numpy as np
import math
from app import config
from app.routing.helper import mymath


def scaleshape(pts, scalingFactor):
    x = pts[:, 0]
    y = pts[:, 1]
    (rho, phi) = mymath.array_cart2pol(x, y)
    (x, y) = mymath.array_pol2cart(rho * scalingFactor, phi)
    return np.c_[x, y]


def scalerings(rings, num):
    for i in range(len(rings)):
        rings[i].bp = rings[i].bp + num
        rings[i].radius = rings[i].bp * config.AXIAL_RISE / (2 * math.pi)
    return rings


# @param: base_layer
#   a list of Ring objects
# @param: num
#   how many layers to add
def reinforce(base_layer, num):
    all_layers = []
    for r in base_layer:
        newRingGroup = RingGroup(r)
        newRingGroup.extend(num)
        all_layers.append(newRingGroup)
    return all_layers


def numberRingGroup(groups):
    num = 1
    for i, g in enumerate(groups, start=1):
        g.groupnum = i
        for r in g.rings:
            r.ringnum = num
            num += 1
    return


def renumberActiveRings(groups):
    num = 1
    for i, g in enumerate(groups, start=1):
        g.groupnum = i
        for r in g.rings:
            if r.active == True:
                r.ringnum = num
                num += 1


def find(planes, ringnum):
    for p in planes:
        for r in p.rings:
            if r.ringnum == ringnum and r.active:
                return r
    raise IndexError("Did not find that ring.")


class RingGroupAngle:
    def __init__(self):
        self.value = 0

    def __get__(self, instance, owner):
        # if the RingGroup is already set
        if len(instance.rings) > 1:
            ring1 = instance.rings[0]
            ring2 = instance.rings[1]

            r1 = ring1.radius
            h1 = ring1.height

            r2 = ring2.radius
            h2 = ring2.height

            (x, y) = np.array([r2 - r1, h2 - h1])
            (rho, phi) = mymath.cart2pol(x, y)
            return mymath.roundnearest(math.degrees(phi), 5)
        # otherwise return the default angle
        else:
            return math.degrees(self.value)

    def __set__(self, instance, value):
        self.value = value
        # only act if there are rings to adjust
        if len(instance.rings) > 1:
            instance.adjust(value)
        # don't do anything if there is only a single layer
        else:
            pass


class RingGroup:
    angle = RingGroupAngle()

    def __init__(self, *rings, groupnum=-1):
        self.rings = [r for r in rings]
        self.groupnum = groupnum

    def extend(self, num):
        # extend deletes any existing rings and makes a new ringgroup with additional rings
        self.radangle = math.radians(self.angle)
        self.rings = [self.rings[0]]
        baseRing = self.rings[0]
        baseRing.refit()
        for n in range(num):
            baseRingCoords = np.array([baseRing.radius, baseRing.height])
            v = np.array([math.cos(self.radangle), math.sin(self.radangle)])
            nextRingCoords = baseRingCoords + v * config.INTERHELICAL * (n + 1)
            newradius = nextRingCoords[0]
            newbp = mymath.roundnearest(mymath.radius2bp(newradius), base=baseRing.numxovers)
            newtbx = baseRing.tbx
            newnumxovers = baseRing.numxovers
            newheight = nextRingCoords[1]
            newRing = Ring(newradius, newtbx, newnumxovers, newheight)
            newRing.resetRing(newbp, newtbx, newnumxovers, newheight)
            newRing.refit()
            self.rings.append(newRing)

    def adjust(self, angle):
        self.radangle = math.radians(angle)
        baseRing = self.rings[0]
        baseRingCoords = np.array([baseRing.radius, baseRing.height])
        v = np.array([math.cos(self.radangle), math.sin(self.radangle)])
        for n, r in enumerate(self.rings[1:]):
            nextRingCoords = baseRingCoords + v * config.INTERHELICAL * (n + 1)
            newradius = nextRingCoords[0]
            newbp = mymath.roundnearest(mymath.radius2bp(newradius), base=baseRing.numxovers)
            newtbx = baseRing.tbx
            newnumxovers = baseRing.numxovers
            newheight = nextRingCoords[1]
            r.resetRing(newbp, newtbx, newnumxovers, newheight)
            r.refit()

    def __str__(self):
        return str([str(r) + "<br>" for r in self.rings])


class Ring(object):
    def __init__(self, bp, tbx, numxovers, height, active=True, ringnum=-1, dirBit=None):
        # log.debug(__name__,
        #           "Creating RING (bp={}, tbx={}, xovers={}, h={}, dir={})".format(bp, tbx, numxovers, height, dirBit))
        # calculate and set properties

        self.tbx = tbx  # turns between crossovers
        self.numxovers = numxovers
        if not config.CIRCUMFERENCE_ROUNDING:
            self.bp = int(bp - (bp % self.numxovers))
        else:
            self.bp = int(bp)
        self.radius = mymath.bp2radius(self.bp)
        self.height = height
        self.bbx = bp / numxovers  # bases between crossovers7
        self.bpt = self.bbx / self.tbx  # bases per turn
        try:
            self.apb = 360.0 / self.bpt  # angle per base
        except:
            raise RuntimeError("bp: {}".format(self.bp))

        self.connected = []  # connectivity list is for mesh processing

        self.ringnum = ringnum  # -1 is not set
        self.active = active
        self.dirBit = dirBit

    def refit(self):
        ## only iterate if BPT is not satisfied
        if not config.MINBPT <= self.bpt <= config.MAXBPT:
            ## Iterate until stable
            tpxrange = np.arange(config.MINTPX, config.MAXTPX + 1, 1)

            stable = False
            while (stable == False):
                # if valid tpx can be applied, ring is stable
                # otherwise, unstable
                for t in tpxrange:
                    bpt = self.bbx / t
                    if bpt <= config.MAXBPT and bpt >= config.MINBPT:
                        self.resetRing(self.bp, t, self.numxovers, self.height)
                        stable = True
                        break
                    else:
                        stable = False
                # if unstable and bpt>12, double the crossovers
                if stable == False and self.bbx / config.MAXTPX > 12:
                    self.resetRing(self.bp, self.tbx, self.numxovers * 2, self.height)
                # if unstable, scale the ring by one unit of numx
                elif stable == False:
                    # print("I'm scaling!")
                    self.resetRing(self.bp + self.numxovers, self.tbx, self.numxovers, self.height)
        else:
            pass

    def setDir(self, dirBit):
        self.dirBit = dirBit

    def resetRing(self, bp, tbx, numxovers, height):
        # print("attempting tbx={}".format(tbx))
        # calculate and set properties
        self.radius = mymath.bp2radius(bp)
        self.tbx = tbx  # turns between crossovers
        self.numxovers = numxovers
        self.bp = bp
        self.height = height
        self.bbx = bp / numxovers  # bases between crossovers
        self.bpt = self.bbx / self.tbx  # bases per turn
        self.apb = 360.0 / self.bpt  # angle per base

    def __str__(self):
        return "{}, {}, {}, {}".format(self.bp, self.tbx, self.numxovers, self.height)

    def data(self):
        return (self.bp, self.tbx, self.numxovers, self.height, self.active, self.ringnum, self.dirBit)

    def setconnected(self, ring):
        self.connected.append(ring)
