from flask import request, make_response, session, render_template
import numpy as np
import itertools

from . import bp
from .session_helpers import require_session_keys
from .utils import edges_to_linkedlist, is_connected
from app.routing.helper.mymath import bp2radius


@bp.route('/connections')
@require_session_keys('circdata', 'ringdata')
def upload_connections():
    circs = session['circdata']
    rings = session['ringdata']
    ringmap = {}
    for i, r in enumerate(rings):
        ringmap[i] = [bp2radius(r[0]), r[1]]
    connections = []
    for ring1, ring2 in itertools.combinations(ringmap, 2):
        pt1 = np.array(ringmap[ring1])
        pt2 = np.array(ringmap[ring2])
        if np.linalg.norm(pt2-pt1) <= 3.0:
            connections.append([ring1, ring2])
    return render_template('upload-connections.html',
                           existing=connections,
                           circleCoords=circs,
                           wizard_step=2)


@bp.route("/upload_connections", methods=["GET", "POST"])
def connections_input():
    if request.method == "POST":
        data = request.json
        graph = edges_to_linkedlist(data)
        if not data:
            resp = make_response("ERROR: Received no input.", 400)
            return resp
        elif not is_connected(graph, len(session['ringdata'])):
            resp = make_response("ERROR: Not all helices are connected.", 400)
            return resp
        else:
            pass
        session["connections"] = data
        return "Continue"
    else:
        raise RuntimeError("This function should only be called when submitting connections data.")
