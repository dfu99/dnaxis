from app import app
from app import config

from flask import render_template, send_file, send_from_directory, session
import os


# loading screen while waiting for process to finish
@app.route('/processing')
def processing():
    return render_template("processing.html")


@app.route('/download/<path:filename>', methods=['GET', 'POST'])
def download_file(filename):
    output_dir = config.JOBS_DIR
    return send_from_directory(output_dir, filename, as_attachment=True, cache_timeout=0)


# gives a download link
@app.route('/results')
def showresults():
    return render_template('results.html')


# download page
@app.route('/download')
def download():
    outputdir = session['wdir']
    try:
        return send_file(os.path.join("..", outputdir, "export.zip"), as_attachment=True, max_age=0,
                     download_name="export.zip")
    except Exception as e:
        print(e)
        jobid = outputdir[5:]
        return render_template("error.html", jobno=jobid)
