"""
This file helps to generate the connections.json and pathway.json for custom shape inputs
"""


# Assumes that this is entered in the correct pathway order
input = [[0, 0, 328, True, 0, 0, 0],
         [0, 2.6, 328, False, 0, 1, 1],
         [0, 5.2, 328, True, 0, 2, 2],
         [0, 7.8, 328, False, 0, 3, 3],
         [0, 10.4, 328, True, 0, 4, 4],
         [0, 13, 328, False, 0, 5, 5],
         [0, 15.6, 328, True, 0, 6, 6],
         [0, 18.2, 328, False, 0, 7, 7],
         [0, 20.8, 328, True, 0, 8, 8],
         [0, 23.4, 328, False, 0, 9, 9],
         [0, 26, 328, True, 0, 10, 10],
         [0, 28.6, 328, False, 0, 11, 11],
         [0, 31.2, 328, True, 0, 12, 12],
         [0, 33.8, 328, False, 0, 13, 13],
         [0, 36.4, 328, True, 0, 14, 14],
         [0, 39, 328, False, 0, 15, 15],
         [0, 41.6, 328, True, 0, 16, 16],
         [0, 44.2, 328, False, 0, 17, 17],
         [0, 46.8, 328, True, 0, 18, 18],
         [0, 49.4, 328, False, 0, 19, 19],
         [0, 52, 328, True, 0, 20, 20],
         [0, 54.6, 328, False, 0, 21, 21],
         ]


'''Make the connections file here'''
import itertools
import numpy as np

# Runs the same automatic connections finder as the UI website
connections = dict([])
connections["connections"] = dict([])

for m1, m2 in itertools.combinations(input, 2):
    pt1 = np.array([m1[0], m1[1]])
    pt2 = np.array([m2[0], m2[1]])
    if np.linalg.norm(pt2 - pt1) <= 2.7:  # Interhelical distance
        try:
            connections["connections"][m1[6]].append(m2[6])
        except:
            connections["connections"][m1[6]] = [m2[6]]
        try:
            connections["connections"][m2[6]].append(m1[6])
        except:
            connections["connections"][m2[6]] = [m1[6]]


# EXPORT CONNECTIONS
import json

cout = dict([])
cout["connections"] = []
for key in connections["connections"]:
    source = [key, 0]
    targs = connections["connections"][key]
    targets = [[x, 0] for x in targs]
    cout["connections"].append({"source": source, "targets": targets})

with open("connections_3.0_0.0.json", 'w') as outfile:
    json.dump(cout, outfile)


'''Make the pathway file here'''
pathway = []

for i in range(38):
    for line in input:
        if line[5] == i:
            pathway.append(line[6])

pout = dict([])
pout["pathway"] = []

for i1, i2 in zip(pathway[:-1], pathway[1:]):
    source = [i1, 0]
    target = [i2, 0]
    pout["pathway"].append({"source": source, "target": target})

with open("pathway.json", 'w') as outfile:
    json.dump(pout, outfile)
