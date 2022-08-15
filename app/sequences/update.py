"""
Created by Daniel Fu (Reif Lab, Duke University) at 10/29/2020

Name        : update
Project     : cadaxisdna
Description : Dynamically allocates scaffold sequences depending on list of active files in directory
                If files change, run refresh() again
Interpreter : Python 3.7.4
"""

import json
import os


def readseq(filename):
    """
    Reads single line sequence from file
    :param filename: string filename
    """
    with open(os.path.join(".", filename), 'r') as f:
        seq = f.readline()
    f.close()
    return seq


def writeseqlib(out_data):
    """
    Writes json-formated dict to .json file
    :param out_data: json-formated dict
    :return: None, generates file
    """
    # Convert to json format
    with open('seqlib.json', 'w') as outfile:
        json.dump(out_data, outfile, separators=(',', ':'))
    outfile.close()
    print("Written data to", 'seqlib.json')


def readseqlib():
    """
    Reads sequence library file
    :return: dict of sequence lengths
    """
    json_data = open(os.path.join('sequences', 'seqlib.json')).read()
    data = json.loads(json_data)
    return data


def refresh():
    """
    Gets a new list of sequences based on what files are available in the folder
    :return: Writes a json with keys=names and values=lengths of scaffolds
    """
    dirname = os.path.join(".")
    jsondata = {}
    for file in os.listdir(dirname):
        if file.endswith(".txt"):
            filename = file.replace(".txt", "")
            jsondata[filename] = len(readseq(file))
    writeseqlib(jsondata)


def sort_seq_lib(lib):
    """
    Sorts the library of sequences by length
    :param lib: dict of sequences
    :return: Sorted library (list by keys)
    """
    sortedlib = sorted(lib.items(), key=lambda x: x["Length"], reverse=True)
    return sortedlib


def seq_remainder(keys, tarlen):
    """
    initiate a dict with scaffold keys each initiated to the target length

    iterate through the list of scaffold sequences

    update each key recursively with scaffold sequences that are shorter than the key itself

    pick the key that has the lowest remainder

    :param keys: dict containing sequence library
    :param tarlen: Target length (int)
    """
    # Initialize each scaffold sequence key with remaining unassigned scaffold
    tardict = {}
    for k in keys:
        tardict[k] = {"Length": tarlen, "Path": [k]}

    # Update, subtract sequence length from each key
    for k in keys:
        tardict[k]["Length"] -= keys[k]

    # Recursively update for each nonnegative key
    for k in tardict.keys():
        if tardict[k]["Length"] > 0:
            tkeys = {}  # Truncate the dict to only recurse on shorter scaffold sequences
            for j in keys:
                if keys[j] < keys[k]:
                    tkeys[j] = keys[j]
            # If remaining is non-empty
            if tkeys:
                try:
                    next_key = seq_remainder(tkeys, tardict[k]["Length"])
                    tardict[k]["Length"] = next_key["Length"]
                    [tardict[k]["Path"].append(p) for p in next_key["Path"]]
                except TypeError:  # No return
                    pass

    # Return minimum remainder (greatest negative value)
    minkey = None
    minlen = min([tardict[t]["Length"] for t in tardict.keys()])
    # Skip positive (invalid) keys
    for k in tardict.keys():
        if minlen <= tardict[k]["Length"] <= 0:
            minkey = k
    try:
        return tardict[minkey]
    except KeyError:
        return None


if __name__ == "__main__":
    refresh()
    # seqlens = readseqlib()
    # print(seqlens)
    # print("The answer is", seq_remainder(seqlens, 10000))
