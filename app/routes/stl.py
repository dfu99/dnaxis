from flask import request, flash, redirect, render_template, url_for, session, current_app
from werkzeug.utils import secure_filename
import os
import datetime
from shapely.geometry import LineString
import numpy as np

from . import bp
from .utils import allowed_file
from app.scaffold_library import SCAFFOLD_LENGTHS
from app.meshing.mycompress import encode, decode
from app.meshing.graph import build_graph, graph_stl
from app.meshing.symmesh import stl2ptcloud, flatten, readstl, fithelices, fitboundary, openshape
from app import config


@bp.route('/stl')
def upload_stl():
    session['fromstl'] = 0
    return render_template('upload-stl.html')


@bp.route('/show-example', methods=["GET", "POST"])
def process_example():
    if request.method == 'POST':
        stl = request.json
        print("request is {}".format(stl))
        print(os.getcwd())
        if stl == "none":
            pass
        elif stl == "gourd":
            stl_filepath = os.path.join('app', 'static', 'stl', 'gourd.stl')
            filename = "gourd"
        elif stl == "vase":
            stl_filepath = os.path.join('app', 'static', 'stl', 'vase.stl')
            filename = "vase"
    return _showstl(stl_filepath, filename)


@bp.route('/upload-stl', methods=["GET", "POST"])
def process_stl():
    if request.method == 'POST':
        if 'file' not in request.files:
            flash('ERROR: No file part')
            return redirect(url_for('main.upload_stl'))
        f = request.files['file']

        if not f.filename == "":
            if not allowed_file(f.filename):
                flash('ERROR: Please check that your file is in Stereolithography (.stl) format')
                return redirect(url_for('main.upload_stl'))
            if f and allowed_file(f.filename):
                wdir = os.path.join(current_app.config['UPLOAD_FOLDER'],
                                    datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
                session['wdir'] = wdir
                if not os.path.exists(wdir):
                    os.makedirs(wdir)
                stl_filepath = os.path.join(wdir, secure_filename(f.filename))
                f.save(stl_filepath)
                filename = f.filename
        else:
            flash('ERROR: No selected file')
            return redirect(url_for('main.upload_stl'))
    else:
        flash('Unexpected Error: Please contact dfu@cs.duke.edu')
        return redirect(url_for('main.upload_stl'))

    return _showstl(stl_filepath, filename)


def _showstl(stl_filepath, filename):
    session['alpha'] = 1
    session['alpha_off'] = 0.5
    session['isOpen'] = False
    session['stl_filepath'] = stl_filepath
    session['filename'] = filename

    mesh = readstl(stl_filepath)
    imgX = graph_stl(mesh, view=(0, 0, 'x'))
    imgY = imgX
    imgZ = imgX
    response_data = {"xview": imgX, "yview": imgY, "zview": imgZ}

    pts = stl2ptcloud(mesh)

    flat_reps = []
    for c in ['x', 'y', 'z']:
        flat_reps.append(flatten(pts, center=c))

    session['flat_reps'] = encode.nparr(flat_reps)

    return response_data


@bp.route('/stl2nodes', methods=['GET', 'POST'])
def stl2nodes():
    if request.method == 'POST':
        opt_axis = 'Z'
        if opt_axis == 'X':
            axis = 0
        elif opt_axis == 'Y':
            axis = 1
        else:
            axis = 2
        flat_reps = decode.nparr(session['flat_reps'])
        pts = flat_reps[axis]

        last_alpha = alpha = session['alpha']
        isOpen = session['isOpen']

        while True:
            try:
                boundary = fitboundary(pts, alpha)
                last_alpha = alpha
            except AttributeError:
                try:
                    boundary = fitboundary(pts, last_alpha)
                    break
                except AttributeError:
                    flash("Sorry! This file could not be parsed correctly.")
                    return redirect((url_for('main.upload_stl')))
            alpha += 1

        active_segment = boundary[0]

        if isOpen:
            active_segment = openshape(active_segment)

        session['halfxsec'] = encode.nparr(np.array(list(active_segment.coords)))
        _ = build_graph(pts=pts, halftrace=active_segment, fulltrace=boundary[1], figsize=(5, 5))

        xsec = LineString(decode.nparr(session['halfxsec']))
        minbp = 72
        mintpx = 2
        numrings = 20
        xovers = 4

        if minbp > mintpx * xovers * config.MINBPT:
            use_tpx = False
        else:
            use_tpx = True

        scafAvailable = 50000

        base_layer = fithelices(user_maxnt=scafAvailable, user_numrings=numrings, user_xoverCount=xovers,
                                        use_tpx=use_tpx, MINCIRCUMFERENCE=minbp, MINTPX=mintpx, xsec=xsec, interhelical=config.INTERHELICAL)

        session['ringdata'] = encode.rings(base_layer)
        print(session['ringdata'])
        session['fromstl'] = 1
        return redirect('/submission')
