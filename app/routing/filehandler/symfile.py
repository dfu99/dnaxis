# -*- coding: utf-8 -*-
"""
Created on Sun Jun  9 20:45:20 2019

@author: Dan
"""

import csv

'''
### IMPORT WORKSPACE FILES
Reads from externally generated (MATLAB) data
Data indicates location and connection of nucleotides
'''


# =============================================================================
# IMPORTS
# =============================================================================
# reads a txt file with symmetrical plane-by-plane origami data
# and converts to float
def getfile(filename):
    data = readfile(filename)
    data = list2DToFloat(data)
    data = format_rings(data)
    return data


def readfile(filename):
    # import innermost base pair ring data for symmetric structures
    # expected csv format: <bases>,<turns between crossovers>,<crossovers>,<height>
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        data = list(reader)
    f.close()
    
    return data


def list2DToFloat(data):
    # convert to float 2-D list
    for i in range(len(data)):
        for j in range(len(data[i])):
            try:
                data[i][j] = round(float(data[i][j]), 6)  # <bases><turns between xovers><crossovers><height>
            except:
                data[i][j] = data[i][j] == 'True'
    return data


def format_rings(data):
    """
    Convert each index to correct format position
    [0] = bp
    [1] = tpx
    [2] = numxovers
    [3] = height
    [4] = blank
    [5] = blank
    [6] = dirbit
    :param data:
    :return:
    """
    newdata = []
    for line in data:
        newdata.append([line[0], line[1], line[2], line[3], -1, -1, line[5]])
    return newdata


# test code
if __name__ == "__main__":
    filename = "../test/Cylinder.csv"
    data = getfile(filename)