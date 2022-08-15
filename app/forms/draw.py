from app import app
from flask import request, session, flash, redirect, url_for
import os
import datetime

from app.routing.helper.mymath import bp2radius
from app import config

from . import scafLen

# re module provides support for regular expressions for checking email address validity
import re

# Make a regular expression for validating an Email
regex = r'\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Z|a-z]{2,}\b'


# Define a function for validating an Email
def check_email(email):
    # pass the regular expression
    # and the string into the fullmatch() method
    if re.fullmatch(regex, email):
        return True
    else:
        return False


# get the STL input file and user input variables
@app.route('/submit', methods=['GET', 'POST'])
def uploader_submission():
    # UPLOAD THE FILE
    if request.method == 'POST':
        # Enforce format requirements
        # Skip email for now
        email = request.form['opt_email']
        if check_email(email):
            session['user-email'] = email
        else:  # email is optional
            session['user-email'] = None
            # flash('Please enter in a valid email address.')
            # return redirect(url_for('submit'))

        # create an instance folder for this job
        # default setting, must always happen
        wdir = os.path.join(app.config['UPLOAD_FOLDER'], datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
        session['wdir'] = wdir
        if not os.path.exists(wdir):
            os.makedirs(wdir)

        # save the input parameters
        session['mintpx'] = request.form['opt_mintpx']
        session['xovercount'] = request.form['opt_xovercount']
        session['scaf'] = request.form['opt_scaffolds']
        if int(request.form.get('opt_lenlow')) >= int(request.form.get('opt_lenup')):
            flash("Invalid staple length bounds.")
            return redirect(url_for('upload_submission'))
        session['lenlow'] = int(request.form.get('opt_lenlow'))
        session['lenup'] = int(request.form.get('opt_lenup'))

        # advanced options
        setattr(config, 'VALIDXOVERTHRESHBP', float(request.form.get('opt_vxotbp')))

        setattr(config, 'VALIDXOVERSPACING_SAME', float(request.form.get('opt_vxoss')))

        setattr(config, 'VALIDXOVERSPACING_ADJ', float(request.form.get('opt_vxosa')))

        # save configurations
        setattr(config, 'MINTPX', int(session['mintpx']))
        setattr(config, 'INTERHELICAL', float(request.form.get('opt_interdist')))
        setattr(config, 'LENLOW', session['lenlow'])
        setattr(config, 'LENUP', session['lenup'])

        # save param string from advanced options
        paramstring = '-debug'

        if request.form.get('opt_shape'):
            paramstring += ' -shape'
        if request.form.get('opt_frs'):
            paramstring += ' -frs'
        if request.form.get('opt_oldrouting'):
            paramstring += ' -oldrouting'
        if request.form.get('opt_savesteps'):
            paramstring += ' -savesteps'
        if request.form.get('opt_fcs'):
            paramstring += ' -fcs'
        if request.form.get('opt_fcm'):
            paramstring += ' -fcm'
        if request.form.get('opt_uvxl'):
            paramstring += ' -uvxl'
        if request.form.get('opt_ox'):
            paramstring += ' -ox'
        if request.form.get('opt_stats'):
            paramstring += ' -stats'
        if request.form.get('opt_valid'):
            paramstring += ' -valid'
        session['paramstring'] = paramstring
        print("Paramstring=", paramstring)

        # save the path to the saved STL file
        session['filename'] = "test"

        # Get the input data
        data = request.form.get('mesh_txt_input')
        data = data.splitlines()
        data = [line.split(",") for line in data]
        data = [[s.replace(" ", "") for s in line] for line in data]
        data = [[int(num[0]), float(num[1]), int(num[2])] for num in data]

        # converts ring objects to plain text to save the ring data
        session['ringdata'] = data

        # error handling: checks if structure surpasses p8064 + phix174
        if "-shape" not in paramstring:
            used_scaf = 0
            for r in session['ringdata']:
                used_scaf += r[0]
            print("used_scaf= ", used_scaf)
            if used_scaf > scafLen['p8064'] + scafLen['phix174']:
                flash('ERROR: That structure will be larger than currently supported scaffold length limits.')
                return redirect(url_for('upload_submission'))

        data = [[bp2radius(int(num[0])), float(num[1])] for num in data]
        session['circdata'] = data

        return redirect('/connections')
    flash("Critical error: Received no data, contact dfu@cs.duke.edu.")
    return redirect(url_for('upload_submission'))