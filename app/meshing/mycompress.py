import io
import zlib
import numpy as np
from . import helixgeom as hxg

class encode:

    def nparr(nparr):
        """
        Returns the given numpy array as compressed bytestring,
        the uncompressed and the compressed byte size.
        """
        bytestream = io.BytesIO()
        np.save(bytestream, nparr)
        uncompressed = bytestream.getvalue()
        compressed = zlib.compress(uncompressed)
        return compressed

    def rings(ringslist):
        data = [r.data() for r in ringslist]
        return data

    def ringGroups(layers):
        planes = []
        for group in layers:
            modules = []
            for ring in group.rings:
                modules.append(ring.data())
            planes.append({'modules':modules, 'groupnum':group.groupnum})
        return planes


class decode:

    def nparr(bytestring):
        """
        Returns the given bytestring as a decompressed numpy array
        """
        return np.load(io.BytesIO(zlib.decompress(bytestring)))

    def rings(datalist):
        rings = [hxg.Ring(int(d[0]), int(d[1]), int(d[2]), float(d[3]), active=d[4], ringnum=d[5], dirBit=d[6]) for d in datalist]
        return rings

    def simple_rings(datalist):
        rings = [hxg.Ring(int(d[0]), 2, 4, int(d[1]), ringnum=i) for i, d in enumerate(datalist)]
        return rings

    def ringGroups(data):
        planes = []
        for line in data:
            modules = line['modules']
            groupnum = line['groupnum']
            rings = [hxg.Ring(m[0], m[1], m[2], m[3], active=m[4], ringnum=m[5], dirBit=m[6]) for m in modules]
            newRingGroup = hxg.RingGroup(*rings, groupnum=groupnum)
            planes.append(newRingGroup)
        return planes