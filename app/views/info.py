from app import app
from flask import render_template, session


# Home page
@app.route('/')
def home():
    session['fromstl'] = 0
    return render_template('index.html')


# Test page
@app.route('/test')
def test():
    return render_template('test.html')


# Documentation page
@app.route('/tutorial')
def docs_page():
    session['fromstl'] = 0
    return render_template('tutorial.html')
