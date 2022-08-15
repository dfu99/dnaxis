from app import app
from flask import request, flash, redirect, url_for, session
from werkzeug.utils import secure_filename
from .utils import allowed_file
import os
import datetime
from shapely.geometry import LineString

from . import scafLen
from app.meshing.mycompress import encode, decode
from app.meshing.graph import build_graph, graph_stl
from app.meshing.symmesh import stl2ptcloud, flatten, readstl, fithelices, fitboundary, openshape
import numpy as np
from app import config

# Receives and processes STL file
@app.route('/upload-stl', methods=["GET", "POST"])
def process_stl():
    # UPLOAD THE FILE
    if request.method == 'POST':
        # check if the post request has the file part
        if 'file' not in request.files:
            flash('ERROR: No file part')
            return redirect(url_for('upload_stl'))
        # read from the uploaded file
        f = request.files['file']

        # check for uploaded file first
        # then check for example file
        # check if it's a valid file
        # if user does not select file, return an error
        if not f.filename == "":  # file was selected
            if not allowed_file(f.filename):  # if file format is incorrect
                flash('ERROR: Please check that your file is in Stereolithography (.stl) format')
                return redirect(url_for('upload1_stl'))
            if f and allowed_file(f.filename):  # if file format is correct and file exists
                # create an instance folder for this job
                # default setting, must always happen
                wdir = os.path.join(app.config['UPLOAD_FOLDER'],
                                    datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
                session['wdir'] = wdir
                if not os.path.exists(wdir):
                    os.makedirs(wdir)
                # save the input STL file
                stl_filepath = os.path.join(wdir, secure_filename(f.filename))
                f.save(stl_filepath)
                filename = f.filename
        else:  # file was NOT selected, flash an error, stay on the page
            flash('ERROR: No selected file')
            return redirect(url_for('upload_stl'))
    else:
        flash('Unexpected Error: Please contact dfu@cs.duke.edu')
        return redirect(url_for('upload_stl'))

    # set a default alpha for boundary fitting later
    session['alpha'] = 1
    # set a default view alpha for mesh multi preview
    session['alpha_off'] = 0.5
    # set a default isOpen for horizontal edges
    session['isOpen'] = False
    # save the path to the saved STL file
    session['stl_filepath'] = stl_filepath
    session['filename'] = filename

    mesh = readstl(stl_filepath)
    imgX = graph_stl(mesh, view=(0, 0, 'x'))
    # TODO: Disabled until rotation is implemented in STL preview
    # imgY = graph_stl(mesh, view=(90, 90, 'y'))
    # imgZ = graph_stl(mesh, view=(90, 0, 'z'))
    imgY = imgX
    imgZ = imgX
    response_data = {"xview": imgX, "yview": imgY, "zview": imgZ}

    pts = stl2ptcloud(mesh)

    # flatten the representation along each axis
    flat_reps = []
    for c in ['x', 'y', 'z']:
        flat_reps.append(flatten(pts, center=c))

    session['flat_reps'] = encode.nparr(flat_reps)

    return response_data


@app.route('/stl2nodes', methods=['GET', 'POST'])
def stl2nodes():
    if request.method == 'POST':
        # Get the selected axis of rotation
        # TODO: Disabled until rotation is implemented in STL preview
        # opt_axis = session['opt_axis'] = request.form['opt_axis']
        opt_axis = 'Z'
        if opt_axis == 'X':
            axis = 0
        elif opt_axis == 'Y':
            axis = 1
        else:
            axis = 2
        # Pick out that particular projection from a list
        flat_reps = decode.nparr(session['flat_reps'])
        pts = flat_reps[axis]

        # Pull settings
        alpha = session['alpha']
        isOpen = session['isOpen']

        # Run the boundary finding algorithm
        while True:
            try:
                boundary = fitboundary(pts, alpha)
                last_alpha = alpha
            except AttributeError:
                boundary = fitboundary(pts, last_alpha)
                break
            alpha += 1

        # [0] = halfxsec
        # [1] = fullxsec
        active_segment = boundary[0]

        # Check whether to fill in horizontal segments
        if isOpen:
            active_segment = openshape(active_segment)

        # Save the unrevolved shape
        session['halfxsec'] = encode.nparr(np.array(list(active_segment.coords)))
        _ = build_graph(pts=pts, halftrace=active_segment, fulltrace=boundary[1], figsize=(5, 5))

        xsec = LineString(decode.nparr(session['halfxsec']))
        # get the input parameters
        # minbp = int(session['minbp'])
        # mintpx = int(session['mintpx'])
        # numrings = int(session['numrings'])
        # xovers = int(session['xovercount'])
        minbp = 72
        mintpx = 2
        numrings = 20
        xovers = 4

        # check if the minimum circumference should overwite the calculated minimum derived from mintpx
        if minbp > mintpx * xovers * config.MINBPT:
            use_tpx = False
        else:
            use_tpx = True

        # determine how much scaffold to use for the fit
        # scaf = session['scaf']
        # customScafLen = request.form.get('opt_scafLen', type=int)
        # if customScafLen != scaf:
        #     scafAvailable = customScafLen
        # else:
        #     scafAvailable = scafLen[scaf]

        # Simply parsing the structure, independent of scaffold usage
        scafAvailable = 50000

        # perform the fit
        base_layer = fithelices(user_maxnt=scafAvailable, user_numrings=numrings, user_xoverCount=xovers,
                                        use_tpx=use_tpx, MINCIRCUMFERENCE=minbp, MINTPX=mintpx, xsec=xsec, interhelical=config.INTERHELICAL)

        # converts ring objects to plain text to save the ring data
        session['ringdata'] = encode.rings(base_layer)
        print(session['ringdata'])
        session['fromstl'] = 1
        return redirect('/submission')
