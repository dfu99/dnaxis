from app import app
from flask import session, render_template


# set scaffold pathway
@app.route('/pathway')
def upload_pathway():
    circs = session['circdata']
    connections = session['connections']
    return render_template('upload-pathway.html', edges=connections, circleCoords=circs)


