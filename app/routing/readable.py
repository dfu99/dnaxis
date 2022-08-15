"""
Created by Daniel Fu (Reif Lab, Duke University) at $(DATE)

Name        : readable.py
Project     : cadaxisdna
Description : Helps convert objects into human readable output
Interpreter : Python 3.7.4
"""

from .helper import strandnav, log


def xovers(xovers):
    """
    Displays list of crossovers between two rings in a readable format
    :param xovers: list of type FullXover
    :return: printout
    """
    if not xovers:
        log.system("\nNo crossovers were listed!")
    for i, x in enumerate(xovers, start=1):
        log.system("\nCrossover {}".format(i))
        (bp, h) = x.n1.up().getposition_round()
        log.system("\t Ring 1 : {:6} {:6} Pos : {:6}".format(bp, h, round(x.n1.pos), 6))
        (bp, h) = x.n2.up().getposition_round()
        log.system("\t Ring 2 : {:6} {:6} Pos : {:6}".format(bp, h, round(x.n2.pos), 6))


def xovers_pap(dnaorigami):
    """
    Print each crossover and whether it is parallel or anti-parallel
    :param dnaorigami:
    :return:
    """
    modules = dnaorigami.get_modules()
    for m in modules:
        for xover in m.objxovers:
            log.system("Anti-parallel crossover: {}"
                  "\n\tn15 is loop: {}"
                  "\n\tn13 is loop: {}"
                  "\n\tn25 is loop: {}"
                  "\n\tn23 is loop: {}".format(xover.is_anti_parallel(),
                  strandnav.checkloop(xover.n1.n5),
                  strandnav.checkloop(xover.n1.n3),
                  strandnav.checkloop(xover.n2.n5),
                  strandnav.checkloop(xover.n2.n3)))


def helixfeatures(helix):
    """
    Displays indices and position of nucleotides on a helix and whether they are a feature
    :param helix: type DNAHelix
    :return: printout
    """
    scaf = helix.scaf
    stap = helix.stap
    for i, (scnucl, stnucl) in enumerate(zip(scaf, stap)):
        log.system(i,
              strandnav.isfeature(scnucl), strandnav.isfeature(stnucl),
              round(scnucl.phi, 4))


def connections(dnaorigami, conns_dict):
    """
    Print out the note for each module to help with making sense of the connections graph
    :param dnaorigami: The Origami object
    :param conns_dict: The connections dictionary
    :return:
    """
    for key in conns_dict.keys():
        source_helix = dnaorigami.planes[key[0]].modules[key[1]].note
        for value in conns_dict[key]:
            target_helix = dnaorigami.planes[value[0]].modules[value[1]].note
            log.system("{}, {}: {} connected to {}, {}: {}.".format(key, round(dnaorigami.planes[key[0]].height, 1),
                                                               source_helix,
                                                               value, round(dnaorigami.planes[value[0]].height, 1),
                                                               target_helix))

    # # Displaying nucl ID and checking for present features on helix
    # helix = dnaorigami.planes[3].modules[0].helix
    # readable.helixfeatures(helix)
    # log.system(helix.up().getsimpleposition())
    #
    # # Displaying crossovers and positions between two rings
    # ring1 = dnaorigami.planes[0].modules[0]
    # ring2 = dnaorigami.planes[1].modules[0]
    #
    # appliedxovers = dnaorigami.get_xovers_btwn_rings(ring1, ring2)
    # ring1xovers = appliedxovers[ring1]
    # ring2xovers = appliedxovers[ring2]
    #
    # readable.xovers(ring1xovers)
    #
    # # Removing crossovers to decouple rings
    # for r in ring1xovers:
    #     r.undo()
    #
    # appliedxovers = dnaorigami.get_xovers_btwn_rings(ring1, ring2)
    # ring1xovers = appliedxovers[ring1]
    # ring2xovers = appliedxovers[ring2]
    #
    # readable.xovers(ring1xovers)

    # # check a specific strand containing specified nucleotide
    # # prints out column
    # log.system("\t# check a specific strand containing specified nucleotide")
    # testnucl = dnaorigami.get_nucl(1201)
    # from routing.helper import strandnav
    # strand = strandnav.getstrand(testnucl)
    # log.system("Strand of nucleotide {}".format(testnucl.numid))
    # for nucl in strand:
    #     log.system(nucl.numid)

    # # check a specific strand containing specified nucleotide
    # # prints out row
    # log.system("\t# check a specific strand containing specified nucleotide")
    # testnucl = dnaorigami.get_nucl(311)
    # from routing.helper import strandnav
    # strand = strandnav.getstrand(testnucl)
    # log.system("Strand of nucleotide {}".format(testnucl.numid))
    # log.system(",".join([str(nucl.numid) for nucl in strand]))

    # # print out adjacent numids
    # for nucl in dnaorigami.planes[2].modules[0].shape.graph:
    #     log.system(nucl.prev_point, nucl.next_point)
    # for nucl in dnaorigami.planes[2].modules[0].helix.stap:
    #     log.system(nucl.numid, nucl.toFive, nucl.toThree)

    # # check every strand
    # all_strands = dnaorigami.get_all_staples()
    # for strand in all_strands:
    #     log.system(",".join([str(nucl.numid) for nucl in strand]))
    #     _ = input("")

    # # prints out row of nucl on each strand
    # for i, p in enumerate(dnaorigami.planes):
    #     for j, m in enumerate(p.modules):
    #         log.system("### dnaorigami.planes[{}].modules[{}].helix.scaf ###".format(i, j))
    #         # log.system(len(m.helix.scaf))
    #         log.system(",".join([str(nucl.numid) for nucl in m.helix.scaf]))
    #         log.system("### dnaorigami.planes[{}].modules[{}].helix.stap ###".format(i, j))
    #         # log.system(len(m.helix.stap))
    #         log.system(",".join([str(nucl.numid) for nucl in m.helix.stap]))

    # # For checking the position of arcs
    # for i, p in enumerate(dnaorigami.planes):
    #     for j, m in enumerate(p.modules):
    #         log.system(i, j, m.shape.center, m.bp)

    # # Check the 5'-3' orientation of each module
    # for i, p in enumerate(dnaorigami.planes):
    #     for j, m in enumerate(p.modules):
    #         nucl = m.helix.scaf[0]
    #         nucl5 = strandnav.goto5(nucl)
    #         nucl3 = strandnav.goto3(nucl)
    #         log.system(i, j, nucl5.center, nucl3.center)

    # For checking angles of crossover nucleotides
    # for p in pathway:
    #     midx1 = p[0]
    #     midx2 = p[1]
    #     module1 = dnaorigami.planes[midx1[0]].modules[midx1[1]]
    #     module2 = dnaorigami.planes[midx2[0]].modules[midx2[1]]
    #     xovers2 = module1.appliedxoversto(module2)
    #
    #     log.system(xovers2[0].n1.coords)
    #     log.system(xovers2[0].n2.coords)
    #
    #     log.system(xovers2[0].n1.n5.theta, xovers2[0].n1.n5.asymtheta)
    #     log.system(xovers2[0].n1.n3.theta, xovers2[0].n1.n3.asymtheta)
    #     log.system(xovers2[0].n2.n5.theta, xovers2[0].n2.n5.asymtheta)
    #     log.system(xovers2[0].n2.n3.theta, xovers2[0].n2.n3.asymtheta)

    # # Check the crossover spacings
    # module1 = dnaorigami.planes[0].modules[4]
    # module2 = dnaorigami.planes[1].modules[4]
    # dnaorigami.edge2edge_crossover_spacing(module1, module2)

    # Compare the desired module to module angle, and the theta and asymtheta of the chosen NuclPairs
    #   seems kind of useless though
    # for key in connections:
    #     module1 = dnaorigami.planes[key[0]].modules[key[1]]
    #     for val in connections[key]:
    #         log.system("Connecting {}, {} to {}, {}".format(key[0], key[1], val[0], val[1]))
    #         module2 = dnaorigami.planes[val[0]].modules[val[1]]
    #         (p1, p2) = makeshape.findclosestpoints(module1.shape, module2.shape)
    #         (x, y, z) = p2.coordinates - p1.coordinates
    #         (rho, theta, phi) = mymath.cart2sph(x, y, z)
    #         phi = math.degrees(phi)
    #         module1.setadjacent(module2, phi)
    #         NuclPairs = module1.getxoversto(module2, 'stap')
    #         for np in NuclPairs:
    #             log.system("{:14} | phi: {:6} | theta: {:6} | asymtheta {:6}".format(
    #                 str(np.numid), phi, np.theta, round(np.asymtheta, 0)))

    # Troubleshooting
    # Plot the node graph in Python to compare with oxViewer
    # import matplotlib.pyplot as plt
    # fig, ax = plt.subplots(figsize=(10, 10))  # note we must use plt.subplots, not plt.subplot
    # for plane in dnaorigami.planes:
    #     for module in plane.modules:
    #         x, y = (np.array([]), np.array([]))
    #         node_graph = module.shape.graph
    #         for point in node_graph:
    #             x = np.append(x, point.coordinates[0])
    #             y = np.append(y, point.coordinates[1])
    #         ax.plot(x, y)
    # dim = 300
    # ax.set_xlim((-dim, dim))
    # ax.set_ylim((-dim, dim))
    # plt.show()


if __name__ == "__main__":
    pass
