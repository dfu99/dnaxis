from app import app
from flask import session, render_template
from app.meshing.mycompress import decode
from app.routing.helper.mymath import bp2radius
from app import config


# upload file
@app.route('/submission')
def upload_submission():
    if session['fromstl']:
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
                           fromstl=session['fromstl'],
                           existing=coords,
                           scaffolds=config.AVAIL_SEQUENCES,
                           interhelical=config.INTERHELICAL,
                           lenlow=config.LENLOW,
                           lenup=config.LENUP,
                           vxotbp=config.VALIDXOVERTHRESHBP,
                           vxoss=config.VALIDXOVERSPACING_SAME,
                           vxosa=config.VALIDXOVERSPACING_ADJ)
