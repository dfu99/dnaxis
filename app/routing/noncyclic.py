"""
Created by Daniel Fu (Reif Lab, Duke University) at 3/12/2021

Name        :
Project     :
Description :
Interpreter : Python 3.7.4
"""


import numpy as np


def is_module_non_cyclic(module):
    """
    Checks whether the module is non-cyclic
    :param module:
    :return:
    """
    for nucl in module.helix.scaf:
        if nucl.toThree == -1 or nucl.toFive == -1:
            return True
    return False


def verify(dnaorigami):
    """
    Checks whether the Origami object has non-cyclic sections
    :param dnaorigami: Origami object
    :return: True/False
    """
    for m in dnaorigami.get_modules():
        if is_module_non_cyclic(m):
            return True
    return False


def gap_exists(m1, m2):
    """
    Checks whether both modules are non-cyclic
    :param m1: modules object
    :param m2: modules object
    :return: True/False
    """
    if is_module_non_cyclic(m1) and is_module_non_cyclic(m2):
        return True
    else:
        return False


def gap_xover_exists(m1, m2):
    """
    Checks whether the modules are connected by a normal scaffold crossover or a gap-spanning crossover
    Does so by checking the separation distance of nucleotides within a NuclPair
    :param m1:
    :param m2:
    :return:
    """
    for xover in m1.objxovers + m2.objxovers:
        n1n5 = xover.n1.n5
        n1n3 = xover.n1.n3
        n2n5 = xover.n2.n5
        n2n3 = xover.n2.n3
        if np.linalg.norm(n1n5.node.coordinates - n1n3.node.coordinates) > 1.0 and \
                np.linalg.norm(n2n5.node.coordinates - n2n3.node.coordinates) > 1.0:
            return True
    return False


if __name__ == "__main__":
    pass