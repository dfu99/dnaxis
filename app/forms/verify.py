from app import app
from flask import request, redirect

# loading screen while waiting for process to finish
@app.route('/verified', methods=["GET", "POST"])
def verified():
    if request.method == 'POST':
        return redirect("/processing")