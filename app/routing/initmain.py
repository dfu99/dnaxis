#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 00:41:44 2018

@author: dfu

initmain.py
main function call deciding between sym or asym functions
split from initorigami.py

CHANGELOG:
    4/2/2018:
        - Add function for importing base layer data of symmetric structures
        - Add function for populating nucleotides of symmetric origami
        - Add function for importing point data of asymmetric structures
        - Add function for population nucleotides of asymmetric origami
    4/18/2018:
        - Consolidate import and draw sym and asym functions here
    4/19/2018:
        - Rename to initorigami.py
        - Functions are for initializing nucleotides to fit imported shape data
    5/11/2018:
        - Add routing classes directly to origami_data in draw_sym()
    7/12/2018:
        - Split import_sym() function further into
            - sym_import - imports ONLY base layer
            - sym_expand - performs any multi-layer calculation
            - data2ring - initializes ring modules according to imported data
        - Split sym and asym code into two different files
        - Rename initorigami.py function to initmain.py
    8/21/2018:
        - Add log output
"""
from .helper import log
from . import initsym


# controller function
# chooses between symmetric or asymmetric motif
def draw_origami(filename, motif):  # DEPRECATED
    foldername = 'input/'
    filename = foldername + filename
    if motif == 's':
        log.out(__name__, "Importing SYMMETRIC routing using RING objects from imported data in", filename+".csv")
        planes = initsym.draw_sym(filename)
    return planes