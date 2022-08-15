from app import app

from flask import session, render_template
import numpy as np
from app.routing.helper.mymath import bp2radius

# set staple connections
@app.route('/connections')
def upload_connections():
    import itertools
    circs = session['circdata']
    rings = session['ringdata']
    ringmap = {}
    for i, r in enumerate(rings):
        ringmap[i] = [bp2radius(r[0]), r[1]]
    connections = []
    for ring1, ring2 in itertools.combinations(ringmap, 2):
        pt1 = np.array(ringmap[ring1])
        pt2 = np.array(ringmap[ring2])
        if np.linalg.norm(pt2-pt1) <=3.0:
            connections.append([ring1, ring2])
    return render_template('upload-connections.html', existing=connections, circleCoords=circs)