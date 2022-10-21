from app import app
from flask import session, render_template
import numpy as np
from app.routing.helper.mymath import bp2radius
import itertools
from app import config

# run some preliminary checks on the submission
@app.route('/pre-submit')
def pre_submit():
    connections = session['connections']
    pathway = session['pathway']
    nodes = session['circdata']
    labels = []
    messages = ""
    numissues = [0 for i in nodes]

    # Check for single layer structures
    message = '<p>Some parts of the structure are single-walled that may have lower stability. ' \
              'If this is intended, please ignore this message.</p>'
    for idx in range(len(nodes)):
        distinct_edges = []
        for edge in connections:
            if idx in edge:
                distinct_edges.append(edge)
        for edge in pathway:
            if idx in edge:
                if edge in distinct_edges or edge[::-1] in distinct_edges:
                    pass
                else:
                    distinct_edges.append(edge)
        if len(distinct_edges) > 2:
            labels.append("#009E73")
        else:
            labels.append("#F0E442")
            numissues[idx] += 1
            if message not in messages:
                messages += "<br>"
                messages += message

    # Check for stacked crossovers
    message = '<p>Some helices may be colliding given the set interhelical distance ({}). ' \
              'While this is sometimes acceptable, errors can more easily occur.</p>'.format(config.INTERHELICAL)
    rings = session['ringdata']
    ringmap = {}
    for i, r in enumerate(rings):
        ringmap[i] = [bp2radius(r[0]), r[1]]
    stacked = []
    iterlist = list((i, j) for ((i, _), (j, _)) in itertools.combinations(enumerate(ringmap.keys()), 2))
    for b in iterlist:
        pt1 = np.array(ringmap[b[0]])
        pt2 = np.array(ringmap[b[1]])
        # Fudge this by about 0.1, because conversion from bps to cartesian is not ever perfect
        if np.linalg.norm(pt2-pt1) < config.INTERHELICAL - 0.1:
            stacked.append(b[0])
            stacked.append(b[1])
    if len(stacked) > 0:
        for i in stacked:
            numissues[i] += 1
        messages += "<br>"
        messages += message

    return render_template('verify.html', circleCoords=nodes, stapEdges=connections, scafEdges=pathway,
                           labels=labels, numissues=numissues, stacked=stacked, messages=messages)