"""
Created by Daniel Fu (Reif Lab, Duke University) at 9/28/2020

Name        : strand.py
Project     : cadaxisdna
Description : Strand functions
Interpreter : Python 3.7.4
"""

import random
from .helper import log


class Mixin:
    def parse_strands(self, modules):
        """
        When adding new modules to the Origami, also add their JoinedStrand objects, ignoring duplicates
        :param modules: modules that are being added
        :return: Adds JoinedStrands to attribute
        """
        if type(modules) == list:
            for m in modules:
                refnucl = m.helix.scaf[0]
                strand = refnucl.get_top_strand()
                self.strands.add(strand)
        else:
            refnucl = modules.helix.scaf[0]
            strand = refnucl.get_top_strand()
            self.strands.add(strand)

    def apply_strand_connections(self, connections):
        """
        Set the adjacency of JoinedStrand objects
        :param connections: Dict
        :return: None
        """
        for key in connections:
            strand1 = self.planes[key[0]].modules[key[1]].helix.scaf[0].get_top_strand()
            for val in connections[key]:
                strand2 = self.planes[val[0]].modules[val[1]].helix.scaf[0].get_top_strand()
                strand1.set_adjacent(strand2)

    def apply_strand_pathway(self, pathway):
        """
        Sets the order of the strands according to the input pathway
        :param pathway: Dict
        :return: None
        """
        edges_list = []
        for edge in pathway:
            mdx1 = edge[0]
            strand1 = self.planes[mdx1[0]].modules[mdx1[1]].helix.scaf[0].get_top_strand()
            if edge[1]:
                mdx2 = edge[1]
                strand2 = self.planes[mdx2[0]].modules[mdx2[1]].helix.scaf[0].get_top_strand()
                edges_list.append((strand1, strand2))
        # log.system("Pathway order in terms of consolidated helices")
        # log.system(edges_list)
        return edges_list

    def shuffle_pathway(self, pathway, connections):
        """
        Changes module pairs for the pathway
        :param pathway: Dict
        :param connections: Dict
        :return: new pathway
        """

        """
        Get the JoinedStrand of the source module
        Get the JoinedStrand of the target module
        Randomly choose another module on the source strand
        Get its adjacent modules from the connections list
        Check which module from the connections list is on the target strand
        That is the new pathway
        """
        newpathway = []
        for edge in pathway:
            mdx1 = edge[0]
            mdx2 = edge[1]
            md1 = self.planes[mdx1[0]].modules[mdx1[1]]
            md2 = self.planes[mdx2[0]].modules[mdx2[1]]
            strand1 = md1.helix.scaf[0].get_top_strand()
            strand2 = md2.helix.scaf[0].get_top_strand()

            _md1 = random.choice(strand1.modules)
            _mdx1 = self.get_module_index(_md1)
            for adj_mdx in connections[_mdx1]:
                adj_m = self.planes[adj_mdx[0]].modules[adj_mdx[1]]
                if adj_m in strand2.modules:
                    _md2 = adj_m
                    _mdx2 = self.get_module_index(_md2)
                    break
            newpathway.append((_mdx1, _mdx2))
        log.system("Shuffling pathway")
        for edge in newpathway:
            mdx1 = edge[0]
            mdx2 = edge[1]
            md1 = self.planes[mdx1[0]].modules[mdx1[1]]
            md2 = self.planes[mdx2[0]].modules[mdx2[1]]
            log.system(md1, mdx1, "/", md2, mdx2)
        return newpathway


if __name__ == "__main__":
    pass