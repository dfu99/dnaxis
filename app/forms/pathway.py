from app import app
from flask import request, make_response, session
from .utils import edges_to_linkedlist, is_connected, has_cycle, edge_len_bound, is_lteq_degree

# Receives pathway data from AJAX
@app.route("/upload_pathway", methods=["GET", "POST"])
def pathway_input():
    if request.method == "POST":
        data = request.json
        # print("data:", data)
        graph = edges_to_linkedlist(data)
        # print("graph:", graph)
        # ERROR: Nothing drawn
        # ERROR: Nodes not connected
        # ERROR: Cycle in graph
        # ERROR: Nodes connected more than degree 2
        # ERROR: Connections <2.2nm, >3.0nm
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