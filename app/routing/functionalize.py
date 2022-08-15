'''

functionalize.py

Post-processing modifications to an origami for adding functional groups at specific positions

'''

from .helper import strandnav, dnaconnector, mymath
from . import noncyclic, makeshape, modules
import numpy as np
from operator import itemgetter


def addBase(strand, base, placement):
    """
    adds the selected base to a strand at the specified placements
    :param strand:
      a list of Nucleotide objects
    :param base:
      the base letter to add. can be any letter
    :param placement:
      a string of placement positions
      '5' = Append to 5' end of the strand
      '3' = Append to 3' end of the strand
      'c' = Append to both endpoints of crossovers
    :return: strand
      returns the strand as string for immediate export
      we do this because we don't want to worry about the 3D placement of modifications
    """
    newstrand = ''
    for n in strand:
        if '5' in placement and n.toFive == -1:
            newstrand += base+n.nucl
        elif '3' in placement and n.toThree == -1:
            newstrand += n.nucl+base
        elif 'c' in placement and strandnav.isXover5End(n):
            newstrand += n.nucl+base
        elif 'c' in placement and strandnav.isXover3End(n):
            newstrand += base+n.nucl
        else:
            newstrand += n.nucl
    return newstrand


def extensions_join(nucl1, nucl2):
    """
    Joins two strands while also filling any spacing with polyT
    :param nucl1:
    :param nucl2:
    :return:
    """
    pass


class Mixin:
    def apply_extensions(self):
        """
        Create manual connections between non-adjacent rings.
        :return: None
        """
        ext_edges = []

        threshold = 4
        # extensions format
        #   dict{key=source modules, value=target modules}
        # For each source module
        for source_index in self.extensions.keys():
            # print("Look for source_index", source_index, type(source_index))
            source_module = self.get_module_by_index(source_index)
            # Get all the source nicks
            source_nicks = source_module.get_nicks()
            # print("Num nicks on {} = {}".format(source_module, len(source_nicks)))

            # For each target module of the source module
            for target_index in self.extensions[source_index]:
                target_module = self.get_module_by_index(target_index)
                # Get all the target nicks
                # Split into 5' and 3' positions
                target_nucl3_all = []
                target_nucl5_all = []
                if noncyclic.is_module_non_cyclic(target_module):
                    for nucl in target_module.helix.stap:
                        # Get the gap endpoints
                        if nucl.__strand3__ == -1:
                            target_nucl3_all.append(nucl)
                        if nucl.__strand5__ == -1:
                            target_nucl5_all.append(nucl)
                else:
                    target_nicks = target_module.get_nicks()
                    for nick in target_nicks:
                        target_nucl3_all.append(nick.nuclpair.n5)
                        target_nucl5_all.append(nick.nuclpair.n3)
                # Calculate distances between source nick and target nick
                # For each source nick
                for sn in source_nicks:
                    sn3 = sn.nuclpair.n3
                    for tn3 in target_nucl3_all:
                        d3 = np.linalg.norm(sn3.center - tn3.center)
                        # If within threshold, save this as a tuple <sn, tn, distance>
                        if d3 < threshold:
                            ext_edges.append((sn3, tn3, d3))
                    sn5 = sn.nuclpair.n5
                    for tn5 in target_nucl5_all:
                        d5 = np.linalg.norm(sn5.center - tn5.center)
                        # If within threshold, save this as a tuple <sn, tn, distance>
                        if d5 < threshold:
                            ext_edges.append((sn5, tn5, d5))
        # Combine all source nicks into one list for checking later
        source_nicks_all = []
        [source_nicks_all.append(s.nuclpair.n3) for s in source_nicks]
        [source_nicks_all.append(s.nuclpair.n5) for s in source_nicks]
        # Sort ext_edges by distance
        srt_ext_edges = sorted(ext_edges, key=itemgetter(2), reverse=False)
        # Create a mapping for source to target nicks
        ext_map = {}
        # For each ext_edge
        for edge in srt_ext_edges:
            # Split the tuple
            sn = edge[0]
            tn = edge[1]
            d = edge[2]
            try:  # Check if source nucleotide already has a value
                # If new min distance AND target nucleotide unused
                if d < ext_map[sn]['distance'] and all(tn != ext_map[x]['target'] for x in list(ext_map.keys())):
                    # Add another mapping of target nick and distance to the source nick key
                    ext_map[sn]['target'] = tn
                    ext_map[sn]['distance'] = d
            except KeyError:  # Source nucleotide not mapped yet
                # Check again for tn collision, otherwise do not add the edge
                if all(tn != ext_map[x]['target'] for x in list(ext_map.keys())):
                    ext_map[sn] = {'target': tn, 'distance': d}

            # Check if all dict keys are set
            if all(sn in list(ext_map.keys()) for sn in source_nicks_all):
                break

        # Apply each extension
        for key in ext_map:
            sn = key
            tn = ext_map[key]['target']
            if sn.toThree == -1 and tn.toFive == -1:
                dnaconnector.nuclconnect(sn, tn)
            elif sn.toFive == -1 and tn.toThree == -1:
                dnaconnector.nuclconnect(tn, sn)
            else:
                # print(ext_map)
                # print(sn.toFive, sn, sn.toThree, tn.toFive, tn, tn.toThree, ext_map[key]['distance'])
                raise RuntimeError("This extension is invalid.")

    def replace_bond(self, n5, n3, len):
        """
        Replaces a bond with a uniformly distributed line of nucleotides
        :param n5: 5' Nucleotide
        :param n3: 3' Nucleotide
        :param len: number of additional nucleotides
        :return:
        """
        # Check if it is .center or .node.center
        dirvec = mymath.unitvector(n3.center - n5.center)
        dist = np.linalg.norm(n3.center - n5.center)
        inc = dist/(len + 2)
        newline = makeshape.Line(len, n5.center, dirvec, increment=inc)
        newmodule = modules.LineModule(len, 1, 0, n5.center[2], True)


