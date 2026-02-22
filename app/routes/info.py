from flask import render_template, session
from . import bp


@bp.route('/')
def home():
    session['fromstl'] = 0
    return render_template('index.html')


@bp.route('/test')
def test():
    return render_template('test.html')


@bp.route('/tutorial')
def docs_page():
    session['fromstl'] = 0
    return render_template('tutorial.html')
