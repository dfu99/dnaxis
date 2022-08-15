"""
Created by Daniel Fu (Reif Lab, Duke University) at 4/18/2021

Name        : extensions
Project     : cadaxisdna
Description : for working with extensions json file format
Interpreter : Python 3.7.4
"""


import csv
import json
from ..helper import log

'''
### IMPORT WORKSPACE FILES
Reads from externally generated manual json
Data indicates desired location of extensions
'''


# =============================================================================
# IMPORTS
# =============================================================================
# reads a txt file with symmetrical plane-by-plane origami data
# and converts to float
def getfile(filename):
    data = readfile(filename)
    data = list_removenull(data)
    data = list_2dtoint(data)
    data = format_extensions(data)
    return data


def readfile(filename):
    # import innermost base pair ring data for symmetric structures
    # expected csv format: <bases>,<turns between crossovers>,<crossovers>,<height>
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        data = list(reader)
    f.close()
    return data


def jsondatafromfile(filename):
    """
    Reads json data
    :param filename: string input filename, must be in same path as script
    :return: json data in json-formated dict
    """
    json_data = open(filename).read()
    data = json.loads(json_data)
    # Convert from json to extensions format
    extensions = {}
    for entry in data['extensions']:
        extensions[tuple(entry["source"])] = [tuple(idx) for idx in entry["targets"]]
    return extensions


def datatojson(xtxn_data, filename):
    """
    Writes json-formated dict to .json file
    :param xtxn_data: extensions format
    :param filename: string output filename
    :return: None, generates file
    """
    # Convert to json format
    out_data = {"extensions": []}
    for key in xtxn_data:
        entry = {"source": key, "targets": xtxn_data[key]}
        out_data["extensions"].append(entry)
    with open(filename, 'w') as outfile:
        json.dump(out_data, outfile, separators=(',', ':'))
    outfile.close()
    log.system("Written data to {}".format(filename))


# removes the empty trailing elements in a non-uniform csv
def list_removenull(li):
    for _li in li:
        while True:
            try:
                _li.remove('')
            except:
                break
    return li


# convert to int 2-D list
def list_2dtoint(data):
    for i in range(len(data)):
        for j in range(len(data[i])):
            data[i][j] = int(data[i][j])  # <adjacent rings>
    return data


def format_extensions(data):
    """
    Formats a 2D list of extensions data to the proper format
    Proper format is Dict
    Keys: numrings (index 0 of list)
    Values: Extended numrings (indices 1: of list)
    :param data: 2D list
    :return: Formatted data structure
    """
    newdata = {}
    for line in data:
        newdata[line[0]] = line[1:]
    return newdata


# test code
if __name__ == "__main__":
    filename = "../test/Cylinderc.csv"
    data = getfile(filename)