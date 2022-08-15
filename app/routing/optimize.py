# -*- coding: utf-8 -*-
"""
Created on Sat May 25 22:35:10 2019

@author: Dan
"""

from app import config
from .helper import dnaconnector, log, strandnav
from .score import seeding


class Mixin:
    # returns a list of strands that are still potential targets for
    # optimization
    def valid_opt_strands(self):
        valid = []
        all_staples = self.get_all_staples()
        # go through all staples in the routing
        for strand in all_staples:
            # look at the 5' or 3' adjacent strand
            end5 = strandnav.goto5(strand[0])
            end3 = strandnav.goto3(strand[0])

            # save this region
            end5region = strandnav.getregion(end5)
            end3region = strandnav.getregion(end3)

            adj5nucl = end5.__strand5__
            adj3nucl = end3.__strand3__

            # condition 1:
            # the strand is not seeded
            cond1 = not seeding.has_seed(strand, config.SEEDLEN)

            # condition 2:
            # either 5' or 3' adjacent strand has a non-adjacent seeded region
            [reg5, reg3] = get_adj_regions(strand)
            # condition 2a:
            # 5' adjacent strand has a non-adjacent seeded region and 
            # the adjacent region has nucleotides to spare
            cond2a = strand_seeded_not_region(adj5nucl) and \
                     len(reg5) - (14 - len(end5region)) > config.SEED_BORROW_LIMIT

            # condition 2b:
            # 3' adjacent strand has a non-adjacent seeded region and 
            # the adjacent region has nucleotides to spare
            cond2b = strand_seeded_not_region(adj3nucl) and \
                     len(reg3) - (14 - len(end3region)) > config.SEED_BORROW_LIMIT

            # look for potential optimization cases
            if cond1 and (cond2a or cond2b):
                valid.append(strand)
        return valid

    """
    POST-PROCESS
    Optimizing nicks
    """

    # bug:
    # will sometimes try to perform optimization using the same adjacent strands
    # this is because the algorithm does a single pass
    # but ignores strands that may have been changed as a result of a previous operation
    # update:
    # needs to re-evaluate after every optimization step
    # catch:
    # check how many strands are left, if the number stops changing
    # we've hit a stuck state
    def optimize_seeds(self):
        # while the optimization has not halted
        while True:
            valid = self.valid_opt_strands()
            log.system("Valid left: {}".format(len(valid)))
            try:
                strand = valid.pop(0)
            except IndexError:
                break

            # find the 5' and 3' ends of the strand
            end5 = strandnav.goto5(strand[0])
            end3 = strandnav.goto3(strand[0])

            # save this region
            end5region = strandnav.getregion(end5)
            end3region = strandnav.getregion(end3)

            # look at the 5' or 3' adjacent nucleotide
            adj5nucl = end5.__strand5__
            adj3nucl = end3.__strand3__

            # get the adjacent 5' or 3' region
            [reg5, reg3] = get_adj_regions(strand)
            # condition 2a:
            # 5' adjacent strand has a non-adjacent seeded region and 
            # the adjacent region has nucleotides to spare
            cond2a = strand_seeded_not_region(adj5nucl) and \
                     len(reg5) - (14 - len(end5region)) > config.SEED_BORROW_LIMIT

            # condition 2b:
            # 3' adjacent strand has a non-adjacent seeded region and 
            # the adjacent region has nucleotides to spare
            cond2b = strand_seeded_not_region(adj3nucl) and \
                     len(reg3) - (14 - len(end3region)) > config.SEED_BORROW_LIMIT

            if cond2a:
                shiftnick(reg5, end5region, 14 - len(end5region))
            elif cond2b:
                shiftnick(reg3, end3region, 14 - len(end3region))
            else:
                pass


def shiftnick(from_reg, to_reg, num):
    """
    Shifts the nick between from_reg and to_reg num nucleotides towards from_reg
    Assumptions:
        * from_reg and to_reg are both 5'-3' strands (ordered list of nucleotides)
        * from_reg and to_reg are adjacent
    :param from_reg: Substrand, List of Nucleotide(object), Source region to move nucleotides
    :param to_reg: Substrand, List of Nucleotide(object), Target region to move nucleotides
    :param num: Number of nucleotides to move
    :return: None
    """
    log.system("Moving {} nucleotides from region {} to region {}".format(num, len(from_reg), len(to_reg)))
    # Check 5'-3' relative directionality of adjacent regions
    if to_reg[0] == from_reg[-1].__strand3__:  # from_reg is 5' adjacent to to_reg
        nick = [from_reg[-1], to_reg[0]]
        new_nick = [from_reg[-1 - num], from_reg[-1 - (num - 1)]]
    elif to_reg[-1] == from_reg[0].__strand5__:  # from_reg is 3' adjacent to to_reg
        nick = [to_reg[-1], from_reg[0]]
        new_nick = [from_reg[0 + num], from_reg[0 + (num - 1)]]
    else:
        raise RuntimeError("There is no matching adjacent pair of nucleotides.")
    # join old nick
    dnaconnector.nuclconnect(nick[0], nick[1])
    # make new nick
    dnaconnector.nuclbreak(new_nick[0])


# returns adjacent regions
def get_adj_regions(strand):
    # look at the 5' or 3' adjacent strand
    end5 = strandnav.goto5(strand[0])
    end3 = strandnav.goto3(strand[0])

    # get the 5' and 3' adjacent nucleotides to the strand
    adj5nucl = end5.__strand5__
    adj3nucl = end3.__strand3__

    # get the nearest adjacent region on the adjacent strand
    reg5 = strandnav.getregion(adj5nucl)
    reg3 = strandnav.getregion(adj3nucl)
    return (reg5, reg3)


def strand_seeded_not_region(nucl):
    """
    Checks if input region has another seeded region
    :param nucl: Region specified by nucleotide in that region
    """
    bOut = False  # initiate output as having not found another seed region

    # get the region denoted by the nucleotide
    ref_region = strandnav.getregion(nucl)
    # get the entire strand
    strand = strandnav.getstrand(nucl)
    # get all regions on the strand
    all_regions = strandnav.get_all_regions(strand)

    # filter for seeded regions that are not the original region
    for region in all_regions:
        if reg_is_seed(region, config.SEEDLEN) and region != ref_region:
            # another region is a seed
            bOut = True
        else:
            # other region is not a seed, does not affect output
            pass
    return bOut


def reg_is_seed(region, seedlen):
    """
    Checks if a region is a seed
    :param region: A list of sequential nucleotides representing a substrand (region) of a longer strand
    :param seedlen: seed length, default 14
    """
    if len(region) >= seedlen:
        return True
    else:
        return False
