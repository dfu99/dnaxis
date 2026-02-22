from flask import request, session, render_template, flash, redirect, url_for, current_app
import os
import datetime
import re

from . import bp
from .session_helpers import require_session_keys
from app.meshing.mycompress import decode
from app.routing.helper.mymath import bp2radius
from app import config
from app.sequences.update import update_custom_seq
from app.scaffold_library import SCAFFOLD_LENGTHS

# Email validation regex
_email_regex = r'\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Z|a-z]{2,}\b'


def _check_email(email):
    if re.fullmatch(_email_regex, email):
        return True
    return False


@bp.route('/submission')
def upload_submission():
    if session.get('fromstl'):
        rings = decode.rings(session['ringdata'])
        coords = []
        dirBit = 1
        for r in rings:
            newline = [bp2radius(r.bp), r.height, dirBit]
            coords.append(newline)
            dirBit = int(not dirBit)
    else:
        coords = []
    return render_template('upload-submit.html',
                           fromstl=session.get('fromstl', 0),
                           existing=coords,
                           scaffolds=config.AVAIL_SEQUENCES,
                           interhelical=config.INTERHELICAL,
                           lenlow=config.LENLOW,
                           lenup=config.LENUP,
                           vxotbp=config.VALIDXOVERTHRESHBP,
                           vxoss=config.VALIDXOVERSPACING_SAME,
                           vxosa=config.VALIDXOVERSPACING_ADJ,
                           wizard_step=1)


@bp.route('/submit', methods=['GET', 'POST'])
def uploader_submission():
    if request.method == 'POST':
        email = request.form['opt_email']
        if _check_email(email):
            session['user-email'] = email
        else:
            session['user-email'] = None

        wdir = os.path.join(current_app.config['UPLOAD_FOLDER'],
                            datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
        session['wdir'] = wdir
        if not os.path.exists(wdir):
            os.makedirs(wdir)

        session['mintpx'] = request.form['opt_mintpx']
        session['xovercount'] = request.form['opt_xovercount']
        if int(request.form.get('opt_lenlow')) >= int(request.form.get('opt_lenup')):
            flash("Invalid staple length bounds.")
            return redirect(url_for('main.upload_submission'))
        session['lenlow'] = int(request.form.get('opt_lenlow'))
        session['lenup'] = int(request.form.get('opt_lenup'))

        setattr(config, 'VALIDXOVERTHRESHBP', float(request.form.get('opt_vxotbp')))
        setattr(config, 'VALIDXOVERSPACING_SAME', float(request.form.get('opt_vxoss')))
        setattr(config, 'VALIDXOVERSPACING_ADJ', float(request.form.get('opt_vxosa')))
        session['VALIDXOVERTHRESHBP'] = float(request.form.get('opt_vxotbp'))
        session["VALIDXOVERSPACING_SAME"] = float(request.form.get('opt_vxoss'))
        session["VALIDXOVERSPACING_ADJ"] = float(request.form.get('opt_vxosa'))

        setattr(config, 'MINTPX', int(session['mintpx']))
        setattr(config, 'INTERHELICAL', float(request.form.get('opt_interdist')))
        setattr(config, 'LENLOW', session['lenlow'])
        setattr(config, 'LENUP', session['lenup'])

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

        session['filename'] = "test"

        data = request.form.get('mesh_txt_input')
        data = data.splitlines()
        data = [line.split(",") for line in data]
        data = [[s.replace(" ", "") for s in line] for line in data]
        data = [[int(num[0]), float(num[1]), int(num[2])] for num in data]

        print("[INFO]:", session["scaf"])
        if 'custom' in session["scaf"]:
            custom_scaf_data = request.form.get('custom_scaf_txt')
            if not update_custom_seq(custom_scaf_data):
                flash("Invalid sequence")
                return redirect(url_for('main.upload_submission'))

        session['ringdata'] = data

        if "-shape" not in paramstring:
            used_scaf = 0
            for r in session['ringdata']:
                used_scaf += r[0]
            print("used_scaf= ", used_scaf)
            if 'custom' in session["scaf"] and used_scaf > len(custom_scaf_data):
                flash('ERROR: The custom scaffold length is not long enough for the designed structure.')
                return redirect(url_for('main.upload_submission'))
            elif used_scaf > SCAFFOLD_LENGTHS['p8064'] + SCAFFOLD_LENGTHS['phix174']:
                flash('ERROR: That structure will be larger than currently supported scaffold length limits.')
                return redirect(url_for('main.upload_submission'))

        data = [[bp2radius(int(num[0])), float(num[1])] for num in data]
        session['circdata'] = data

        return redirect('/connections')
    flash("Critical error: Received no data, contact dfu@cs.duke.edu.")
    return redirect(url_for('main.upload_submission'))


@bp.route("/update_sequences", methods=["GET", "POST"])
def update_sequences():
    if request.method == "POST":
        data = request.json
        session['scaf'] = data
        return "Continue"
    else:
        raise RuntimeError("Error in receiving refreshing use_sequences.")
