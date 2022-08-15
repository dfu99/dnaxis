#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  8 20:03:06 2018

@author: dfu

export.py

CHANGELOG:
    - 9/8/2018:
        - File created. Separated export functions from origami.py
"""
import os
from . import functionalize

from .helper import log, strandnav, mymath


# Printing helper function
# Helps to avoid class typing error between int and nucleotide
def adjbase(adj):
    if adj == -1:
        return -1
    else:
        return adj.printid


class Mixin:
    # =============================================================================
    # Output
    # =============================================================================

    def calc_boxsize(self):
        maxbp = 0
        for plane in self.planes:
            for m in plane.modules:
                maxbp = max(m.bp, maxbp)
        return maxbp

    def save_input(self):
        """
        Save the ring list input for reference
        :param filename:
        :return: writes file to disk
        """
        input_filename = os.path.join(self.outputdir, 'modules.csv')
        os.makedirs(os.path.dirname(input_filename), exist_ok=True)
        with open(input_filename, 'w') as f:
            for p in self.planes:
                for r in p.modules:
                    f.write('{:d},{:d},{:f},{:f},{},{}\n'.format(
                        int(r.bp), int(r.turns), r.height, r.bpt, r.dirBit, self.get_module_index(r)))
        f.close()
        log.out(__name__, "Input rings exporting to CSV format.")

    def save_conn(self):
        """
        Save the connections for reference
        :param filename:
        :return: writes file to disk
        """
        conn_filename = os.path.join(self.outputdir, 'connections.csv')
        os.makedirs(os.path.dirname(conn_filename), exist_ok=True)
        with open(conn_filename, 'w') as f:
            for key in self.man_connections:
                f.write('{}'.format(key))
                numvalues = len(self.man_connections[key])
                for n in range(numvalues):
                    f.write(',{}'.format(self.man_connections[key][n]))
                f.write('\n')

        f.close()
        log.out(__name__, "Connections exported to CSV format.")

    def save_path(self):
        """
        Save the pathway for reference
        :param filename:
        :return: writes file to disk
        """
        path_filename = os.path.join(self.outputdir, 'pathway.csv')
        os.makedirs(os.path.dirname(path_filename), exist_ok=True)
        with open(path_filename, 'w') as f:
            for value in self.man_path:
                f.write('{}\n'.format(value))
        f.close()
        log.out(__name__, "Pathway exported to CSV format.")

    def save_asym_path(self):
        """
        Save the pathway for reference. Specific to asym sample files.
        :param filename:
        :return: writes file to disk
        """
        path_filename = os.path.join(self.outpudir, 'pathway.csv')
        os.makedirs(os.path.dirname(path_filename), exist_ok=True)
        with open(path_filename, 'w') as f:
            for value in self.man_path:
                f.write('{}\n'.format(value))
        f.close()
        log.out(__name__, "Pathway exported to CSV format.")

    def print_csv(self):
        csv_filename = os.path.join(self.outputdir, 'sequences.csv')
        os.makedirs(os.path.dirname(csv_filename), exist_ok=True)
        all_strands = self.get_all_strands()

        with open(csv_filename, 'w') as f:
            for s in all_strands:
                end5 = strandnav.goto5(s[0])
                end3 = strandnav.goto3(s[0])
                for n in s:  # string of nucleotides
                    f.write('{0}'.format(n.nucl))
                f.write(',{0}'.format(len(s)))  # length of strand
                # f.write(',{0},{1}'.format(end5.up().up().getsimpleposition(),
                #                           end3.up().up().getsimpleposition()))  # 5' and 3' endpoint rings
                f.write('\n')
        f.close()
        log.out(__name__, "Strand sequences exported to CSV format.")

    def print_functionalized_csv(self, base, options):
        fcsv_filename = os.path.join(self.outputdir, 'funcseq.csv')
        os.makedirs(os.path.dirname(fcsv_filename), exist_ok=True)
        all_strands = self.get_all_strands()
        with open(fcsv_filename, 'w') as f:
            for s in all_strands:
                strand = functionalize.addBase(s, base, options)
                f.write('{},{}\n'.format(strand, len(strand)))
        f.close()
        log.out(__name__, "Functionalized strand sequences export to CSV format.")

    def apply_print_id_to_strand(self, initID, strand):
        for printID, nucl in enumerate(strand, start=initID):
            nucl.printid = printID

    def apply_all_print_id_to_strand(self, allstrands):
        printIDcounter = 0
        for strand in allstrands:
            self.apply_print_id_to_strand(printIDcounter, strand)
            strandlen = len(strand)
            printIDcounter += strandlen

    def print_oxDNA(self, filename="prova", print_mesh=True, print_crossovers=False):
        """
        Generates a oxDNA format output of the origami
        :return: None. Saves output to disk.
        """

        # set filenames and check the disk location
        conf_filename = os.path.join(self.outputdir, filename + '.conf')
        os.makedirs(os.path.dirname(conf_filename), exist_ok=True)
        top_filename = os.path.join(self.outputdir, filename + '.top')
        os.makedirs(os.path.dirname(top_filename), exist_ok=True)
        # xoverd_filename = 'output/'+filename+'/'+'xoverd'+'.xdat'
        # os.makedirs(os.path.dirname(xoverd_filename), exist_ok=True)
        if print_mesh:
            mesh_filename = os.path.join(self.outputdir, filename + '.mesh')
            os.makedirs(os.path.dirname(mesh_filename), exist_ok=True)
        if print_crossovers:
            xovers_filename = os.path.join(self.outputdir, filename + '_xovers.json')
            os.makedirs(os.path.dirname(xovers_filename), exist_ok=True)

        # get all the strands
        all_strands = self.get_all_strands()
        self.apply_all_print_id_to_strand(all_strands)
        boxsize = int(self.calc_boxsize() * 1.1)

        baseIDtotal = 0
        for plane in self.planes:
            for module in plane.modules:
                baseIDtotal = baseIDtotal + len(module.helix.scaf) + len(module.helix.stap)
        strandnumtotal = len(all_strands)

        top = open(top_filename, 'w')
        top.write(f'{baseIDtotal:d} {strandnumtotal:d}\n')
        conf = open(conf_filename, 'w')
        conf.write(
            't = 0\nb = {boxX:f} {boxY:f} {boxZ:f}\nE = 0.000000 0.000000 0.000000\n'.format(boxX=boxsize, boxY=boxsize,
                                                                                             boxZ=boxsize))
        # xoverd = open(xoverd_filename,'w')
        if print_mesh:
            mesh = open(mesh_filename, 'w')
        if print_crossovers:
            labels = []

        for S, strand in enumerate(all_strands, start=1):
            for this_base in strand:
                # write the topology
                A = adjbase(this_base.toThree)
                B = adjbase(this_base.toFive)
                top.write(str(S) + ' ' + this_base.nucl + ' ' + str(A) + ' ' + str(B) + '\n')
                # write the configuration
                [p1, p2, p3] = this_base.position
                [b1, b2, b3] = this_base.backbonevec
                [n1, n2, n3] = this_base.normvec
                conf.write(
                    str(p1) + " " + str(p2) + " " + str(p3) + " " +
                    str(b1) + " " + str(b2) + " " + str(b3) + " " +
                    str(n1) + " " + str(n2) + " " + str(n3) +
                    " 0 0 0 0 0 0\n")
                # write the xover distance spread
                # xoverd.write("{}".format(nicking.distfromxover(this_base))+"\n")
                this_ring = this_base.up().up()
                if print_mesh:
                    mesh.write("{} {} {}\n".format(this_ring.bp, this_ring.height, this_base.numid))
                if print_crossovers:
                    if strandnav.isxover(this_base):
                        labels.append(1)
                    else:
                        labels.append(0)
        top.close()
        conf.close()
        if print_mesh:
            mesh.close()
        if print_crossovers:
            import json
            with open(xovers_filename, 'w') as fxovers:
                data = dict([])
                data['RMSF (nm)'] = labels
                json.dump(data, fxovers)

        return

        log.out(__name__, "Topology and configuration generated to oxDNA format.")

    def print_oxDNA_htrap(self, filename):
        htrap_filename = os.path.join(self.outputdir, filename + '/external.conf')
        os.makedirs(os.path.dirname(htrap_filename), exist_ok=True)

        all_strands = self.get_all_strands()

        htrap = open(htrap_filename, 'w')
        for strand in all_strands:
            for nucl in strand:
                [c1, c2, c3] = nucl.center
                [d1, d2, d3] = nucl.backbonevec
                htrap.write(
                    "{\n" +
                    "type = trap\n" +
                    "particle = {}\n".format(nucl.printid) +
                    "pos0 = {}, {}, {}\n".format(c1, c2, c3) +
                    "stiff = 1.0\n" +
                    "rate = 0.\n" +
                    "dir = {},{},{}\n".format(d1, d2, d3) +
                    "}\n")
        log.out(__name__, "oxDNA harmonic traps file generated.")

    def print_caDNAno(self, filename):
        cadnano_filename = os.path.join(self.outputdir, filename + '.json')
        os.makedirs(os.path.dirname(cadnano_filename), exist_ok=True)

        # initiate the dict
        data = {}
        data["name"] = filename + '.json'
        data["vstrands"] = []

        size = max([m.bp for m in self.get_modules()])
        size = size + (32 - size % 32)
        # give each nucleotide a cadnano lattice position
        idx_pathway = mymath.edgelist2nodelist(self.pathway)
        node_pathway = [self.get_module_by_index(idx) for idx in idx_pathway]
        for num, r in enumerate(node_pathway, start=0):
            r.set_cadnano(size, num)

        #        for i,y in enumerate(self.grid, start=1):
        #            for j,x in enumerate(y, start=23):
        #                x.set_cadnano(size,num)
        #                num+=1

        # set the dict json values
        for i, y in enumerate(self.grid, start=2):
            for j, x in enumerate(y, start=24):
                helix = {}  # clear the dict in each iteration
                helix['row'] = i
                helix['col'] = j
                helix['num'] = x.cdna_num
                helix['scafLoop'] = []
                helix['stapLoop'] = []
                helix['skip'] = [0 for i in range(size)]
                helix['loop'] = [0 for i in range(size)]
                castdata = x.format_cadnano(size)
                helix['scaf'] = castdata['scaf']
                helix['stap'] = castdata['stap']
                helix['stap_colors'] = castdata['colors']
                data["vstrands"].append(helix)
        data['vstrands'] = sorted(data['vstrands'], key=lambda kv: kv['num'])
        DataToJson(data, cadnano_filename)

        return data


def DataToJson(write_data, filename):
    import json
    with open(filename, 'w') as outfile:
        json.dump(write_data, outfile, separators=(',', ':'))
    outfile.close()
    log.system("Written data to {}.".format(filename))
