# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 19:58:19 2019

@author: Dan
pathfile.py

import a manual pathway file
"""

import csv
import json
from ..helper import log

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
    data = list_2dtoint(data)
    data = [x for line in data for x in line]
    return data


def readfile(filename):
    # import innermost base pair ring data for symmetric structures
    # expected csv format: <bases>,<turns between crossovers>,<crossovers>,<height>
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        data = list(reader)
    f.close()
    return data


# convert to int 2-D list
def list_2dtoint(data):
    for i in range(len(data)):
        for j in range(len(data[i])):
            data[i][j] = int(data[i][j])  # <adjacent rings>
    return data


def instructionsfromjson(filename):
    """
    Read pathway routing instructions from json
    :param filename: input file
    :return: instructions data formatted dict
    """
    json_data = open(filename).read()
    data = json.loads(json_data)
    instructions = []
    for edge in data["instructions"]:
        source = edge["source"]
        source["Strand"] = tuple(source["Strand"])
        target = edge["target"]
        target["Strand"] = tuple(target["Strand"])
        entry = (edge["source"], edge["target"])
        instructions.append(entry)
    return instructions


def jsondatafromfile(filename):
    """
    Reads json data
    :param filename: string input filename, must be in same path as script
    :return: json data in json-formated dict
    """
    json_data = open(filename).read()
    data = json.loads(json_data)
    # Convert from json to pathway format
    pathway = []
    for entry in data["pathway"]:
        pathway.append((tuple(entry["source"]), tuple(entry["target"])))
    return pathway


def datatojson(path_data, filename):
    """
    Writes json-formated dict to .json file
    :param path_data: json-formated dict
    :param filename: string output filename
    :return: None, generates file
    """
    # Convert to json format
    out_data = {"pathway": []}
    for edge in path_data:
        entry = {"source": edge[0], "target": edge[1]}
        out_data["pathway"].append(entry)
    with open(filename, 'w') as outfile:
        json.dump(out_data, outfile, separators=(',', ':'))
    outfile.close()
    log.system("Written data to {}.".format(filename))


# test code
if __name__ == "__main__":
    filename = "../test/Cylinderp.csv"
    data = getfile(filename)
