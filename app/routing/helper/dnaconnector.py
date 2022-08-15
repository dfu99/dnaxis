# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 14:43:36 2018

@author: Dan

connect.py

basic manipulation functions for editing 3D DNA files

definitions of
    scafconnect
    stapconnect
    breakstrand
    breakscaf
    
Changelog:
    3/30/2018:
        - Fix directionality. switch each 'cw' with 'ccw' and vice versa
        - Change break value from 0 to -1
    4/1/2018:
        - Change 'CW' and 'CCW' to True, False
    4/19/2018:
        - Adjust for new modular format
    5/13/2018:
        - Modify for new modular code
    7/17/2018:
        - breakstrand and breakscaf are obsolete
        - replaced with nuclbreak
    7/27/2018:
        - deleted obsolete breakstrand, breakscaf, scafconnect and stapconnect functions
    8/17/2018:
        - Added log output
    9/29/2018:
        - Rename nuclconnect as xoverconnect (possibly obsolete)
"""
from . import log

'''
### CONNECTOR FUNCTIONS
'''

# function nuclconnect
# Forms a crossover connection between pairs of nucleotide pairs
# Input two nucleotides (each is the 5' base) of two adjacent strands
# Assumes nucl1 and nucl2 are diagonal positions on adjacent antiparallel helices
# 
# ======xo=======
# ======ox=======
# ======xo=======
# OR:
# ======53=======
# ======35=======
# ======53=======
# 
# where x's are nucl1 or nucl2
def xoverconnect(nucl1,nucl2):
    log.debug(__name__,'Made a crossover at',nucl1.numid,',',nucl2.numid)
    nucl1_5 = nucl1
    nucl1_3 = nucl1.toThree
    nucl2_5 = nucl2
    nucl2_3 = nucl2.toThree
    
    nucl1_5.toThree = nucl2_3
    nucl1_3.toFive = nucl2_5
    nucl2_5.toThree = nucl1_3
    nucl2_3.toFive = nucl1_5
    return


def nuclconnect(nucl5, nucl3):
    """
    Connects two strands at the given adjacent nucleotides
    :param nucl5: 5' Nucleotide(object)
    :param nucl3: 3' Nucleotide(object)
    :return: None
    """
    log.debug(__name__, "Connecting", nucl5.numid, "to", nucl3.numid, ".")
    nucl5.toThree = nucl3
    nucl3.toFive = nucl5


def nuclbreak(nucl):
    """
    Breaks a strand at the given 5' nucleotide position from the 3' nucleotide
    :param nucl: 5' Nucleotide(object)
    :return: None
    """
    log.debug(__name__, 'Broke strand at', nucl.toThree.toFive.numid, ',', nucl.toThree.numid)
    nucl.toThree.toFive = -1
    nucl.toThree = -1
    return


def trueconnect(nucl5, nucl3):
    """
    Connects two co-plane modules
    :param nucl5: 5' Nucleotide(object)
    :param nucl3: 3' Nucleotide(object)
    :return: None
    """
    log.debug(__name__, "Absolute connecting", nucl5.numid, "to", nucl3.numid, ".")
    nucl5.__strand3__ = nucl3
    nucl3.__strand5__ = nucl5


def truebreak(nucl):
    """
    Breaks two co-plane modules, usually if one is deleted
    :param nucl: 5' Nucleotide(object)
    :return: None
    """
    log.debug(__name__, 'Absolute connections broken between', nucl.numid, ',', nucl.__strand3__.numid)
    nucl.__strand3__.__strand5__ = -1
    nucl.__strand3__ = -1
