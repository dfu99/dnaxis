from app import app
from flask import request, make_response, session
from .utils import edges_to_linkedlist, is_connected

# Receives connections data from AJAX
@app.route("/upload_connections", methods=["GET", "POST"])
def connections_input():
    if request.method == "POST":
        data = request.json
        # print("data:", data)
        graph = edges_to_linkedlist(data)
        # print("graph:", graph)
        # WARNING: Connections >3.0 nm
        # ERROR: Nodes not connected
        if not data:
            resp = make_response("ERROR: Received no input.", 400)
            return resp
        elif not is_connected(graph, len(session['ringdata'])):
            resp = make_response("ERROR: Not all helices are connected.", 400)
            return resp
        else:
            pass
        # Otherwise, save and continue
        session["connections"] = data
        return "Continue"
    else:
        raise RuntimeError("This function should only be called when submitting connections data.")
