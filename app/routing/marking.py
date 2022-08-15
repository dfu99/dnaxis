#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  9 19:18:23 2018

@author: dfu

marking.py

Changelog:
    9/9/2018:
        File created
"""
from .helper import log, strandnav
from app import config
from . import modules


class Mixin:
    # =============================================================================
    # Marking functions
    # =============================================================================
    # function init_protected
    # resets the protected list
    def init_protected(self):
        log.out(__name__, "Protected list initialized.")
        self.protected = set([])

    # function visited
    # adds visited nucleotides into the visited list
    def add_protected(self, nucleotide):
        try:
            self.protected.add(nucleotide)
        except:
            log.out(__name__, "ERROR: No visited list found. Initialize it first.")

    # TODO: This only works for crossovers. Extend to work with nicks. Replace nicking.protect_nicks
    def protect(self, nucl, num):
        """
        function protect
        protects num nucleotides around the input 5' nucleotide object, inclusive of input nucl
        :param nucl:
        :param num:
        :return:
        """
        # start with the (5', 3') pair
        nucl3 = nucl.__strand3__
        nucl5 = nucl
        try:
            log.debug(__name__, 'Protecting {} bases around {}-{}.'.format(num, nucl5.numid, nucl3.numid))
        except AttributeError:
            pass
        self.add_protected(nucl3)
        self.add_protected(nucl5)
        self.add_protected(nucl3.Comp)
        self.add_protected(nucl5.Comp)
        # move outward num-1 spaces, such that the total number of nucleotides
        # visited from each 5' or 3' nucleotide outward is equal to num
        for i in range(num - 1):
            nucl3 = nucl3.__strand3__
            if type(nucl3) != int:
                self.add_protected(nucl3)
                self.add_protected(nucl3.Comp)
            else:
                break
        for i in range(num - 1):
            nucl5 = nucl5.__strand5__
            if type(nucl5) != int:
                self.add_protected(nucl5)
                self.add_protected(nucl5.Comp)
            else:
                break

    def reprotect(self):
        """
        Redo the protected list if crossovers have been edited
        :return:
        """
        self.init_protected()
        for nucl in self.get_all_helix_nucl():
            if strandnav.isxover(nucl):
                self.protect(nucl, config.PROTECT)

    def label_module_indices(self):
        """
        Refresh labels on modules if inputs are custom and make_bundle_path was skipped
        :return:
        """
        all_modules = self.get_modules()
        for m in all_modules:
            if type(m) == modules.RingModule or type(m) == modules.ArcModule:
                index = self.get_module_index(m)
                m.note = "({}, {}) {}".format(index[0], index[1], "NA")
