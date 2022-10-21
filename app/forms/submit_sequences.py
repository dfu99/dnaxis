from app import app
from app import config
from flask import request, make_response, session


# Receives use_sequences from AJAX, draw page
@app.route("/update_sequences", methods=["GET", "POST"])
def update_sequences():
    if request.method == "POST":
        data = request.json
        session['scaf'] = data
        print(data)
        print(config.AVAIL_SEQUENCES)
        resp = make_response("Notice: Received data {}, type:{}.".format(data, repr(type(data))), 400)
        return resp
        return "Continue"
    else:
        raise RuntimeError("Error in receiving refreshing use_sequences.")
