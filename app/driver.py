"""
Created by Daniel Fu (Reif Lab, Duke University) at May 13, 2022

Name        : driver.py
Project     : cadaxisdna
Description : Generalized process for creating any shape, usually asymmetric shapes
Interpreter : Python 3.7.4
"""

from .routing.helper import dnaconnector, log
from .routing.filehandler import connfile, pathfile, xtxnfile, symfile
from .sequences import update as sequpd
from .routing import noncyclic, modules, sym_origami
from app import config
from app.routing.sequence import setscafsdir, setscafs
import numpy as np
import os

# adjust some of the settings as needed
# this is to be done before the origami package is initiated because some of the settings are set on import
setscafsdir(os.path.join('app', 'sequences'))

global JOBS_DIR
JOBS_DIR = os.path.join("test")

global TEST_DIR
TEST_DIR = os.path.join("app", "static", "blanks")


def blank_origami(name, output_dir):
    """
    Manually sets up an empty Origami object
    :return: Created Origami object
    """
    # first we import the rings, pathway, and connection data
    sfile = os.path.join(TEST_DIR, 'blank.csv')
    pfile = os.path.join(TEST_DIR, 'blankp.csv')
    cfile = os.path.join(TEST_DIR, 'blankc.csv')
    rings = symfile.getfile(sfile)
    path = pathfile.getfile(pfile)
    connections = connfile.getfile(cfile)
    # then create the origami according to the inputs
    # origami requires rings, pathway, connections
    # this initiates the object with the settings but has not populated the object nucleotides
    dnaorigami = sym_origami.Origami(rings, path, connections, name, output_dir)
    return dnaorigami


def connections_file_exists(filename):
    if os.path.exists(filename):
        return True
    return False


def pathway_file_exists(filename):
    if os.path.exists(filename):
        return True
    return False


def extension_file_exists(filename):
    if os.path.exists(filename):
        return True
    return False


def driver(filename, output_dir, shape, crossover_factor, connect_options, xover_mod,
           auto_scaf_options=False,
           shape_only=False,
           skip_routing=False,
           skip_nicks=False,
           skip_routing_scaf=False,
           skip_routing_stap=False,
           skip_sequence=False,
           force_shuffle_pathway=False,
           force_rand_seq=False,
           force_reseeding=False,
           optimize_xovers=False,
           use_extensions=False,  # Anchors
           add_uvxlinking=False,  # Thymine
           force_xover_density=False,
           replace_long_bonds=False,
           enable_validate=False,
           enable_stats=False,
           enable_debugger=False,
           force_debug_procedures=False,
           force_clean_merge=False,
           force_clear_seam=False,
           use_old_routing=False,
           twist_normalized=False,
           save_steps=False
           ):
    """
    *** Directory parameters
    :param filename: filename prefix
    :param output_dir: directory
    *** Shape and helix parameters
    :param shape: Shape object defined in custom.py
    :param size: helix bundle square cross-section dimensions
    :param crossover_factor: base integer multiple for crossovers
    :param forced_edits: overcompensation for arcs
    :param connect_options: center to center helix distance and percentage of overlap required to be considered nearest
        neighbor helices
    :param xover_mod: additional crossovers to place per nearest neighbor

    **** Driver options
    :param auto_scaf_options: automatically picks scaffolds and number of nicks
    :param shape_only: only creates the helices
    :param skip_routing: skips scaffold, staple routing
    :param skip_routing_scaf: only skips scaffold routing
    :param skip_routing_stap: only skips staple routing
    :param skip_nicks: skips scaffold, staple nicking
    :param skip_sequence: skips sequence application
    :param force_shuffle_pathway: moves pathway along the JoinedStrand (only relevant for asymmetric)
    :param force_rand_seq: random sequence applied to scaffold strands
    :param force_reseeding: greedy algorithm to try and make more seed (default >14 nt) regions for staple binding
    :param optimize_xovers: heuristic algorithm for increasing crossover density and lowering adjacent crossover spacing
    :param use_extensions: reads input file to force connections or handle strands between modules
    :param add_uvxlinking: adds thymines into nicks and crossover for uv crosslinking
    :param force_xover_density: enables xover_mod to be sued
    :param replace_long_bonds: replaces bonds above a certain threshold with polyT spacers (not implemented)
    :param enable_validate: runs tests against the structure
    :param enable_stats: shows some numbers
    :param enable_debugger: does nothing right now
    :param force_debug_procedures: enables some hackish post-processing procedures to clean up nicks and strand lengths
    :param force_clean_merge: forcefully merges short strands, ignoring upper bound
    :param force_clear_seam: removes staple and scaffold crossovers that form seams
    :param use_old_routing: reverts to fixed, old staple, scaffold routing code
    :param twist_normalized: adjusts rings to as close to 10.6 as possible
    :param save_steps: save each iteration of simulated annealing crossover steps

    ***** Output
    Exports oxDNA files configuration and topology, sequence csv, module csv, pathway csv, and connections csv
    :return: None
    """

    log.system(str(locals()))
    # Check debugger variable
    if enable_debugger:
        setattr(config, 'DEBUG_MSGS', True)
    # Initialize origami
    dnaorigami = blank_origami(filename, output_dir)

    # Initiate module creation
    dnaorigami.planes = []
    idx = 0

    """___ ADD ALL MODULES ___"""

    setattr(config, 'TWIST_NORMALIZED', twist_normalized)
    setattr(config, 'XOVER_FACTOR', crossover_factor)
    """Draw the shapes"""
    outputs = shape.construct(idx)
    if type(outputs) == tuple:  # Planes and modules simultaneously defined
        log.system("Planes and modules simultaneously defined.")
        new_planes = outputs[0]
        new_modules = outputs[1]
        for m in new_modules:
            dnaorigami.parse_strands(m)
        for p in new_planes:
            dnaorigami.add_plane(p)

    else:  # Hierarchy not set, module was externally defined
        log.system("Creating planes for defined modules.")
        new_modules = outputs
        for m in new_modules:
            dnaorigami.add_module(m)
            dnaorigami.parse_strands(m)

        for p in dnaorigami.planes:
            log.system("Plane created {}".format(p.height))

    scaf_size = dnaorigami.get_scaf_size()
    log.system("Created modules totaling {} bps.".format(scaf_size))

    # Relabel modules if dealing with custom Ring Module input (Klein bottle)
    dnaorigami.label_module_indices()

    # DEBUG:
    # for m in dnaorigami.get_modules():
    #     m.report_breaks()

    """Automatically determine scaffold settings"""
    if auto_scaf_options:
        # Automatically determine nicking option and which sequences to use
        seqlens = sequpd.readseqlib()
        seq_list = sequpd.seq_remainder(seqlens, scaf_size)
        log.system("Automatically setting scaffold sequence to {}".format(seq_list["Path"]))
        setscafs(*seq_list["Path"])
        if len(seq_list["Path"]) == 1:
            log.system("Automatically configuring to single asymmetric nick.")
            setattr(config, 'SCAF_NICKING', 'asymmetric_single')
        elif len(seq_list["Path"]) == 2:
            log.system("Automatically configuring to double asymmetric nick.")
            if config.ROUTING == 'continuous':
                config.SCAF_NICKING = 'asymmetric_single'
            else:
                config.SCAF_NICKING = 'asymmetric_double'
        else:
            log.system("Automatically configuring to auto nick.")
            config.SCAF_NICKING = 'auto'

    """    Manual crossover adjustment (for 'fixed' XOVERMODE only)    """
    if force_xover_density:
        # manually adjust number of crossovers
        for plane in dnaorigami.planes:
            for module in plane.modules:
                if type(module) == modules.ArcModule:
                    module.numxovers = module.numxovers + xover_mod
                elif type(module) == modules.LineModule:
                    module.numxovers = 1
                else:
                    raise ValueError("Not working.")

    """___IMPORT FILE INSTRUCTIONS TO SET MODULE RELATIONSHIPS ___"""
    # Initialize protected nucleotides set
    dnaorigami.init_protected()

    """    Connections    """
    # File handling strings
    root_connfilename = "connections"
    connfilename = root_connfilename
    for o in connect_options:
        connfilename = connfilename + "_" + str(o)
    filename = os.path.join(output_dir, connfilename + ".json")

    # Options
    # Distance threshold is acceptable bond length of crossovers
    # Overlap is a cutoff for adjacent helices that do not have 100% interfacial alignment
    distance_threshold = connect_options[0]
    overlap_pct = connect_options[1]

    # Skip if shape only parameter requestioned
    if shape_only:
        connections = {}
        pass
    # Otherwise automatically parse connections if not already known
    elif not connections_file_exists(filename):
        # Checks for the 0.0 file

        default_conn_filename = os.path.join(output_dir, root_connfilename + "_" +
                                             str(distance_threshold) + "_" +
                                             str(0.0) + ".json")

        # Prune from the 0.0 version
        if connections_file_exists(default_conn_filename) and overlap_pct > 0:
            log.system("Found the baseline connections and pruning that for connections data fitted to the input "
                       "overlap% of {}.".format(overlap_pct))
            connections = connfile.jsondatafromfile(default_conn_filename)
            connections = dnaorigami.prune_connections(connections, overlap_pct, dist_thresh=2.7)
            connfile.datatojson(connections, filename)
        # Automatically parse for connections using nearest neighbors algorithms using distance and overlap options
        else:
            log.system("Did not find any existing file containing connections data. "
                       "Automatically finding nearest neighbor helices. "
                       "This make take a while.")
            connections = dnaorigami.auto_connections(distance_threshold)
            if overlap_pct > 0:  # Prune again if desired
                connections = dnaorigami.prune_connections(connections, overlap_pct, dist_thresh=2.7)
            connfile.datatojson(connections, filename)
    # Otherwise use existing connections
    else:
        log.system("Using existing connections data from file.")
        connections = connfile.jsondatafromfile(filename)

    # for i, plane in enumerate(dnaorigami.planes):
    #     for k, module in enumerate(plane.modules):
    #         print(module.note, i, k, module)

    # Set connections attributes and all the edges
    log.system("Setting adjacency list for crossovers per each helix.")
    dnaorigami.man_connections = connections
    dnaorigami.apply_asym_connections(connections)
    # Parse the input connections into adjacency of entire helices containing multiple modules
    dnaorigami.apply_strand_connections(connections)

    # # DEBUG
    # log.system("Connections linked list")
    # for s in dnaorigami.strands:
    #     pass
    #     log.system("Joined Strand {} to {}".format(s, s.adjacent))
    #
    # # DEBUG
    # log.system("Distances per edge:")
    # for m in dnaorigami.get_modules():
    #     for ma in m.adjacent:
    #         print("Estimated vector from {} {:6.2f} to {} {:6.2f}: Short Distance: {} | Avg Distance: {} | Angle: {}".
    #               format(m.note, m.get_top_plane().yaw,
    #                      ma.note, ma.get_top_plane().yaw,
    #                      round(m.shortest_dist_to(ma), 2),
    #                      round(m.avg_dist_to(ma), 2),
    #                      round(m.avg_angle_to(ma), 2)))

    """    Pathway    """
    filename = os.path.join(output_dir, "pathway.json")
    if shape_only:
        pathway = []
    # Use pathway instructions if pathway not already parsed
    elif not pathway_file_exists(filename):
        log.system("Did not find any existing file containing pathway data. "
                   "Using the pathway instructions file.")
        pathway = []
        # Translate manual instructions of text description of each module
        # Format as source strand to target strand
        pathway_instructions = pathfile.instructionsfromjson(os.path.join(output_dir, "pathway_instructions.json"))
        for edge in pathway_instructions:
            node1 = edge[0]
            node2 = edge[1]
            for key, plane in enumerate(dnaorigami.planes):
                for value, module in enumerate(plane.modules):
                    if str(node1["Strand"]) in module.note and \
                            node1["Direction"] in module.note and \
                            round(node1["Height"], 2) == round(plane.height, 2):
                        n1key = key
                        n1val = value
                    elif node1["Strand"] == dnaorigami.get_module_index(module) and \
                            node1["Direction"] == "NA" and \
                            round(node1["Height"], 2) == round(plane.height, 2):
                        n1key = key
                        n1val = value
            for key, plane in enumerate(dnaorigami.planes):
                for value, module in enumerate(plane.modules):
                    if str(node2["Strand"]) in module.note and \
                            node2["Direction"] in module.note and \
                            round(node2["Height"], 2) == round(plane.height, 2):
                        n2key = key
                        n2val = value
                    elif node2["Strand"] == dnaorigami.get_module_index(module) and \
                            node2["Direction"] == "NA" and \
                            round(node2["Height"], 2) == round(plane.height, 2):
                        n2key = key
                        n2val = value
            # log.log_debug("(({}, {}), ({}, {}))".format(n1key, n1val, n2key, n2val))
            pathway.append(((n1key, n1val), (n2key, n2val)))
        pathfile.datatojson(pathway, filename)
    # Othwerise use existing path
    else:
        log.system("Using existing pathway data from file.")
        pathway = pathfile.jsondatafromfile(filename)
    log.system("Setting pathway edges for scaffold crossovers.")
    if force_shuffle_pathway:
        # Shuffle the pathway following the same helices as the input pathway
        pathway = dnaorigami.shuffle_pathway(pathway, connections)

    # XOR connections and pathway graphs to de-couple pathway edges from what has been set for connections
    pathway_connections = {}
    for edge in pathway:
        if edge[1] not in connections[edge[0]]:
            try:
                pathway_connections[edge[0]].append(edge[1])
            except KeyError:
                pathway_connections[edge[0]] = [edge[1]]
            try:  # Again, because connections must be explicitly bi-directional
                pathway_connections[edge[1]].append(edge[0])
            except KeyError:
                pathway_connections[edge[1]] = [edge[0]]

    dnaorigami.apply_asym_connections(pathway_connections)
    dnaorigami.apply_strand_connections(pathway_connections)

    dnaorigami.man_path = pathway
    dnaorigami.pathway = pathway
    dnaorigami.pathway_edges = dnaorigami.apply_asym_pathway_edges(pathway)
    dnaorigami.pathway_strands = dnaorigami.apply_strand_pathway(pathway)
    # Debug
    # log.system("Pathway edges are,")
    # for edge in dnaorigami.pathway_edges:
    #     log.system(edge.directed['from'], edge.directed['to'])
    # print("Pathway nodes are,")
    # print(pathway)

    """    Extensions    """
    filename = os.path.join(output_dir, "extensions.json")
    if shape_only or not use_extensions:
        extensions = []
    # Check that file exists
    elif not extension_file_exists(filename):
        raise FileNotFoundError("Extensions are enabled but no instructions file was found.")
    # Otherwise use existing connections
    else:
        log.system("Getting extensions data from file.")
        extensions = xtxnfile.jsondatafromfile(filename)
    dnaorigami.extensions = extensions
    # Debug
    # log.system("Extensions are:")
    # print(extensions)
    # for value in dnaorigami.extensions:
    #     log.system(str(value))

    """___ ROUTING ___"""
    if not shape_only:
        # Test for non-cyclic structures
        log.system("There are non-cyclic segments: {}".format(noncyclic.verify(dnaorigami)))

        """Crossover processing"""
        if not skip_routing:
            if not use_old_routing:
                if not skip_routing_scaf:
                    log.system("Routing scaffold.")
                    dnaorigami.route_asym_scaffold()  # Core step
                    # Clean up any excess bases truncated because of gaps
                    for m in dnaorigami.get_modules():
                        m.clean()
                if not skip_routing_stap:
                    log.system("Routing staples.")
                    dnaorigami.route_asym_staples()  # Core step
                    log.system("Score of the proposed crossover set: {}".format(dnaorigami.xover_set.score()))
            else:  # Legacy methods that do not run heuristics but instead use largest min spacing algorithm
                if not skip_routing_scaf:
                    log.system("Routing scaffold (Legacy).")
                    dnaorigami.route_scaffold()
                if not skip_routing_stap:
                    log.system("Routing staples (Legacy).")
                    dnaorigami.route_staples()


        """Optimization heuristic"""
        # Take the crossover stats and see where additional crossovers can be added
        if optimize_xovers:
            log.system("Applying crossover optimization heuristic algorithm.")
            # Crossover heuristic that adds crossovers into large regions
            log.system("Initial score: {}".format(dnaorigami.xover_set.score()))
            dnaorigami.optimize_crossovers(save_steps=save_steps)
            log.system("Final score: {}".format(dnaorigami.xover_set.score()))

        # Checks for long bonds
        if replace_long_bonds:
            long_bonds = []
            all_xovers = dnaorigami.get_all_xovers()
            for xover in all_xovers:
                n1c = xover.n1.coords
                n2c = xover.n2.coords
                if np.linalg.norm(n2c - n1c) > config.VALIDXOVERTHRESHBP + 5 * config.AXIAL_RISE:
                    long_bonds.append(xover)

            log.system("DEBUG: Num long bonds: {}".format(len(long_bonds)))

        """Nicking processing"""
        if not skip_nicks:
            log.system("Applying nicks.")
            dnaorigami.nick_all()  # Core step
            # TODO: Get rid of this gap hardcoding, this is HORRIBLE coding style
            for gnp in dnaorigami.gapnuclpairs:
                dnaconnector.nuclbreak(gnp.n3.Comp)

        if use_extensions:
            log.system("Setting extended strands.")
            dnaorigami.apply_extensions()

        """Scaffold sequence application"""
        if not skip_sequence:
            log.system("Setting scaffold sequence.")
            if not force_rand_seq:
                dnaorigami.apply_scaf_seq()  # Core step
            else:
                dnaorigami.force_random_seq()  # Core step

        """Post-processing heuristics"""
        if force_reseeding:
            log.system("Optimizing seeds.")
            dnaorigami.optimize_seeds()

        if force_debug_procedures:
            log.system("[WARNING] Running forced manual procedures!")

            # Checks junctions
            # dnaorigami.junction_xovers()  # Only stats, not yet fully implemented

            # First tries to rebalance short strands with adjacent neighbors
            # dnaorigami.rebalance_short_strands()

            if force_clean_merge:
                # Then joins any remaining strands
                dnaorigami.clean_merge()

            if force_clear_seam:
                # Cleans seams
                dnaorigami.clear_seam()

        """ASSESSMENT OF SHAPE INTEGRITY"""
        if enable_validate:
            log.system("Running validity checks.")
            dnaorigami.validate()

            if not dnaorigami.validate_staple_high_bound() or not dnaorigami.validate_staple_low_bound():
                raise RuntimeError("Critical error. Restart.")

    """
    ___ STATS PRINTOUT SUMMARY OF STRUCTURE ___
    """
    if enable_stats:
        pass
        ###
        # Prints the number of crossovers
        log.system("Number of crossovers = {}".format(dnaorigami.num_xovers()))

        ###
        # Prints the number of crossovers as saved on modules
        log.system("Number of saved crossovers = {}".format(len(dnaorigami.get_all_xovers())))

        ###
        # Prints the number of scaffold crossovers
        log.system("Number of scaffold crossovers = {}".format(dnaorigami.num_scaf_xovers()))

        ###
        # Prints the number of staple crossovers
        log.system("Number of staple crossovers = {}".format(dnaorigami.num_stap_xovers()))

        ###
        # For every edge, prints the number of crossovers
        for jstrand1 in dnaorigami.strands:
            for jstrand2 in jstrand1.adjacent:
                log.system("Number of crossovers between {} and {} = {}".format(
                    jstrand1, jstrand2,
                    dnaorigami.num_xovers_by_jstrand(jstrand1, jstrand2)
                ))

        ###
        # For every edge, prints the crossover spacings
        real_score = 0
        for jstrand1 in dnaorigami.strands:
            for jstrand2 in jstrand1.adjacent:
                edge_score = dnaorigami.crossover_spacing_by_jstrand(jstrand1, jstrand2)
                real_score = max(max(edge_score), real_score)
                log.system("Crossover spacings between {} and {}: {}".format(
                    jstrand1, jstrand2, edge_score))
        log.system("Final score: {}".format(real_score))

        ###
        # Prints the number of modules
        num_modules = 0
        for plane in dnaorigami.planes:
            for m in plane.modules:
                num_modules += 1
        log.system("Number of modules = {}".format(num_modules))

        ###
        # Prints which crossovers are parallel or anti-parallel
        # readable.xovers_pap(dnaorigami)

        ###
        # Prints the JoinedStrands
        log.system("Number of helices = {}".format(len(dnaorigami.strands)))

        ###
        # Prints short strands
        # dnaorigami.debug_short_strands()

        ###
        # Outputs to file all crossover bond lengths
        dnaorigami.print_long_bonds(dnaorigami, output_dir)

    # Add thymines into sequences if UV crosslinking is enabled
    if add_uvxlinking:
        pass  # not yet implemented
    # Save the output
    log.system("Exporting output")
    dnaorigami.save()

    return dnaorigami