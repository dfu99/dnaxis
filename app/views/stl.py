from app import app
from flask import session, render_template




# optional tools: start from STL
@app.route('/stl')
def upload_stl():
    session['fromstl'] = 0
    return render_template('upload-stl.html')