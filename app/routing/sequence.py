#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 23:15:32 2018

@author: dfu

CHANGELOG:
    - 9/6/2018: File created

"""
from . import nucleotide
from .helper import log, strandnav
from app import config
import os

# =============================================================================
# File I/O sub-functions
# =============================================================================


# ROUTING CONFIGURATION METHODS
def setscafs(*seq):
    """Changes the order of scaffold sequences that are applied to a structure
    :param seq: new sequence of scaffolds. does not have to include all scaffolds
    """
    current_sequence_order = config.USESEQUENCES
    newsequences = []
    for s in seq:
        newsequences.append(s)
        current_sequence_order.remove(s)
    for s in current_sequence_order:
        newsequences.append(s)
    config.USESEQUENCES = newsequences
    # print("Sequence order set to {}".format(str(newsequences)))


def setscafsdir(dir):
    """Updates the directory to find the sequences
    :param dir: path to sequences/"""
    config.SEQUENCESDIR = dir


# reads single line sequence from file
def readseq(filename):
    with open(os.path.join(config.SEQUENCESDIR, filename+'.txt'), 'r') as f:
        seq = f.readline()
    f.close()
    return seq


# getseq controls for various scaffold sequence options
def getseq(seqname='m13mp18'):
    try:
        return readseq(seqname)
    except FileNotFoundError:
        log.out(__name__,"ERROR: No such sequence. Defaulting to m13mp18.")
        return readseq('m13mp18')


def chooseseqs(size):
    """
    Picks as many scaffold strands as are necessary to complete the routing
    TODO: pick it smartly in order to minimize unused scaffold lengths
    :param size: Number of remaining nucleotides
    :return: List of applied sequences
    """
    chosen = []
    seqnames = iter(config.USESEQUENCES)
    while size > 0:
        try:
            seqname = next(seqnames)
        except StopIteration:
            raise RuntimeError("ERROR: There is not enough scaffold to create the routing. "
                               "No strand output can be generated. Aborting.")
        chosen.append(seqname)
        size -= len(getseq(seqname))
    return chosen


def setseq(init_nucl, seqname):
    """
    Assigns a scaffold sequence 'seq' beginning at 'initnucl'
    :param init_nucl: Nucleotide to begin sequence
    :param seqname: Sequence choice
    :return:
    """

    seq = getseq(seqname)
    log.out(__name__, "Applying {} base pair sequence of length {} to "
                      "scaffold strand starting at {} of length {}.".format(
        seqname, len(seq), init_nucl.numid, strandnav.lenstrand(init_nucl)))
    this_nucl = init_nucl
    cnt = 0
    # For each letter in the sequence 'seq', traverse and apply to the strand 5'-3' starting at this_nucl
    for n in seq:
        this_nucl.nucl = n
        this_nucl.Comp.nucl = nucleotide.comp_nucl(n)
        cnt += 1
        this_nucl = this_nucl.toThree
        if this_nucl == -1:
            log.out(__name__, 'Successfully applied {} base pairs.'.format(cnt))
            return  # Exits after reaching the end of the strand
    # Error handler if somehow the partitioning of crossovers does not match the sequence length
    raise RuntimeError("Scaffold {}, length {} was used up before reaching 3' end.".format(seqname, len(seq)))


class Mixin:
    # =============================================================================
    #  Top level apply scaffold sequence function
    # =============================================================================
    def apply_scaf_seq(self):
        """
        For each start, apply scaffold sequences in the order specified from the configuration
        :return: None
        """
        starts = self.get_scaf_starts()
        log.out(__name__, "Applying {} scaffold sequences. Run apply_scaf_seq to apply different scaffold.".format(
            len(starts)
        ))
        for start, seqname in zip(starts, config.USESEQUENCES):
            setseq(start, seqname)
        
    # =============================================================================
    # Applying scaffold sequence sub-functions
    # =============================================================================
    def get_scaf_starts(self):
        """
        Get all the start positions on the scaffold
        :return: List of all 5' endpoints of the scaffold strands
        """
        # look for a 5' nick
        starts = []
        for module in self.get_modules():
            for nucl in module.helix.scaf:
                if nucl.toFive == -1:
                    starts.append(nucl)
        # ERRORS
        if not starts:  # Scaffold was not nicked and starts is empty
            log.log_error("Scaffold does not have a defined start position. Maybe this was run out of order.")
        num_partitions = len(self.partition_from_removed(self.removed))
        num_starts = len(starts)
        if num_starts != num_partitions:
            log.system("Different number of scaffold start positions ({}) than origami segments ({})!".format(
                num_starts, num_partitions))
            log.log_warning("Wrong number of scaffold nicks.")
        return starts

    def force_random_seq(self):
        """
        Makes sure that all double helices in the structure are complementary but with a random sequence
        :return: None
        """
        for plane in self.planes:
            for module in plane.modules:
                for scnucl in module.helix.scaf:
                    newbase = nucleotide.randnucl()
                    scnucl.nucl = newbase
                    compbase = nucleotide.comp_nucl(newbase)
                    scnucl.Comp.nucl = compbase
