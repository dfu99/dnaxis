from flask import request, make_response, session, render_template

from . import bp
from .session_helpers import require_session_keys
from .utils import edges_to_linkedlist, is_connected, has_cycle, edge_len_bound, is_lteq_degree


@bp.route('/pathway')
@require_session_keys('circdata', 'connections')
def upload_pathway():
    circs = session['circdata']
    connections = session['connections']
    return render_template('upload-pathway.html',
                           edges=connections,
                           circleCoords=circs,
                           wizard_step=3)


@bp.route("/upload_pathway", methods=["GET", "POST"])
def pathway_input():
    if request.method == "POST":
        data = request.json
        graph = edges_to_linkedlist(data)
        if not data:
            resp = make_response("ERROR: Received no input.", 400)
            return resp
        elif not is_connected(graph, len(session['ringdata'])):
            resp = make_response("ERROR: Not all helices are connected.", 400)
            return resp
        elif not is_lteq_degree(graph, 2):
            resp = make_response("ERROR: A helix in the diagram has too many connections.", 400)
            return resp
        elif has_cycle(graph):
            resp = make_response("ERROR: Scaffold routing should not have a cycle.", 400)
            return resp
        elif not edge_len_bound(data, 2.2, 3.0):
            resp = make_response("ERROR: There are helices that are linked by at potentially unstable distances "
                                 "(2.6+/-0.4 nm).", 400)
            return resp
        else:
            pass
        session["pathway"] = data
        return "Continue"
    else:
        raise RuntimeError("This function should only be called when submitting connections data.")
