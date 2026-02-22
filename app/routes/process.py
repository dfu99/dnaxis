from flask import session, render_template, current_app
from flask_mail import Message

from app import config
from app.routing.sequence import setscafs
from app.routing.helper import log, mymath
from app.shapeutils import custom
from app.driver import driver

import os
import zipfile
import datetime
import json
import traceback

from . import bp


def _get_string_options(s, p):
    if p in s:
        return True
    else:
        return False


def _compress_results():
    outputdir = session['wdir']
    try:
        output_file = os.path.join(outputdir, "export.zip")
        if os.path.exists(output_file):
            os.remove(output_file)
        files = ['prova.top', 'prova.conf', 'sequences.csv', 'modules.csv']
        with zipfile.ZipFile(output_file, 'w') as zipf:
            for f in files:
                zipf.write(os.path.join(outputdir, f), f)
        zipf.close()
    except:
        raise RuntimeError


def _save_input(outputdir, ringdata):
    with open(os.path.join(outputdir, "input.txt"), 'w') as f:
        for line in ringdata:
            bps = line[0]
            height = line[1]
            dirbit = line[2]
            f.write("{},{},{}\n".format(bps, height, dirbit))
    f.close()


def _email_error():
    mail = current_app.extensions['mail']
    output_dir = os.path.join(session['wdir'])
    jobid = output_dir[5:]
    if session['user-email']:
        msg = Message('[DNAxiS] job failed', sender='dnaxis.webmaster@gmail.com', recipients=[session['user-email']])
        msg.body = "Rendering of your DNA nanostructure encountered a problem.\n" \
                   "Please contact dfu@cs.duke.edu to submit a bug or if you need assistance." \
                   "Your job number was {}.".format(jobid)
        mail.send(msg)
    return render_template("error.html", jobno=jobid)


def _email_success():
    mail = current_app.extensions['mail']
    output_dir = os.path.join(session['wdir'])
    jobid = output_dir[5:]
    if session['user-email']:
        download_url = config.DOWNLOAD_URL + "{}/export.zip".format(jobid)
        msg = Message('[DNAxiS] job completed', sender='dnaxis.webmaster@gmail.com', recipients=[session['user-email']])
        msg.body = "Rendering of your DNA nanostructure was completed.\n" \
                   "Results can be downloaded from the following link.\n\n" \
                   "{}\n\n" \
                   "Thank you for using DNAxiS!".format(download_url)
        mail.send(msg)


@bp.route('/process')
def upload_process():
    filename = 'test'

    f = open(os.path.join(session['wdir'], "generator_log"), "w")
    pad = config.AXIAL_RISE
    output_dir = os.path.join(session['wdir'])
    f.close()

    log.new('blank', output_dir,
            console=config.LOG_CONSOLE, debug=config.LOG_DEBUG,
            developermode=config.LOG_DEV, log=config.LOG_LOG)
    print("\n")
    log.system("Output to : {}".format(output_dir))
    job_started = datetime.datetime.now()

    connections = session['connections']
    pathway = session['pathway']

    connmap = {}
    for line in connections:
        try:
            connmap[(line[0], 0)].append([line[1], 0])
        except KeyError:
            connmap[(line[0], 0)] = []
            connmap[(line[0], 0)].append([line[1], 0])
        try:
            connmap[(line[1], 0)].append([line[0], 0])
        except KeyError:
            connmap[(line[1], 0)] = []
            connmap[(line[1], 0)].append([line[0], 0])

    connjson = {"connections": []}
    for key in connmap.keys():
        connjson["connections"].append({"source": list(key), "targets": connmap[key]})

    fconnjson = open(os.path.join(output_dir, 'connections_3.0_0.0.json'), 'w')
    json.dump(connjson, fconnjson)
    fconnjson.close()

    pathjson = {"pathway": []}
    for line in pathway:
        pathjson["pathway"].append({"source": [line[0], 0], "target": [line[1], 0]})

    fpathjson = open(os.path.join(output_dir, 'pathway.json'), 'w')
    json.dump(pathjson, fpathjson)
    fpathjson.close()

    rings = session['ringdata']
    shape_class = custom.CustomInput(rings, filename)

    connections_opt = (3.0, 0.0)
    crossover_factor = config.XOVER_FACTOR
    offset_crossover_density = 1
    opt_console_setseq = session['scaf']
    opt_console_scaf_nicking = 'auto'
    paramstring = session['paramstring']

    setscafs(*opt_console_setseq)
    config.SCAF_NICKING = opt_console_scaf_nicking

    config.VALIDXOVERTHRESHBP = session['VALIDXOVERTHRESHBP']
    config.VALIDXOVERSPACING_SAME = session["VALIDXOVERSPACING_SAME"]
    config.VALIDXOVERSPACING_ADJ = session["VALIDXOVERSPACING_ADJ"]

    try:
        driver(filename,
               output_dir,
               shape_class,
               crossover_factor,
               connections_opt,
               offset_crossover_density,
               auto_scaf_options=_get_string_options(paramstring, '-aso'),
               shape_only=_get_string_options(paramstring, '-shape'),
               skip_routing=_get_string_options(paramstring, '-skroute'),
               skip_nicks=_get_string_options(paramstring, '-sknicks'),
               skip_sequence=_get_string_options(paramstring, '-skseq'),
               force_shuffle_pathway=_get_string_options(paramstring, '-fsp'),
               force_rand_seq=_get_string_options(paramstring, '-frs'),
               force_reseeding=_get_string_options(paramstring, '-freseed'),
               optimize_xovers=_get_string_options(paramstring, '-ox'),
               use_extensions=_get_string_options(paramstring, '-ext'),
               add_uvxlinking=_get_string_options(paramstring, '-uvxl'),
               force_xover_density=_get_string_options(paramstring, '-fxd'),
               replace_long_bonds=_get_string_options(paramstring, '-longbonds'),
               enable_validate=_get_string_options(paramstring, '-valid'),
               enable_stats=_get_string_options(paramstring, '-stats'),
               enable_debugger=_get_string_options(paramstring, '-debug'),
               force_debug_procedures=_get_string_options(paramstring, '-debug'),
               force_clear_seam=_get_string_options(paramstring, '-fcs'),
               force_clean_merge=_get_string_options(paramstring, '-fcm'),
               use_old_routing=_get_string_options(paramstring, '-oldrouting'),
               twist_normalized=_get_string_options(paramstring, '-tn'),
               save_steps=_get_string_options(paramstring, '-savesteps'))
        try:
            _compress_results()
        except RuntimeError:
            log.system("Could not collect all the files for export.zip.")
            f.close()
            raise RuntimeError
    except Exception as e:
        _email_error()
        with open(os.path.join(output_dir, "error.txt"), 'a') as f:
            f.write(str(e))
            f.write(traceback.format_exc())
        log.system("Encountered an unknown error. If user supplied an email, they were notified.")
        f.close()
    finally:
        _save_input(output_dir, rings)
        g = open(os.path.join(output_dir, "settings"), "w")
        g.write("[SHAPE]\n")
        g.write("NAME={}\n".format(filename))
        for param in shape_class.inputs:
            g.write("{}={}\n".format(param, shape_class.inputs[param]))
        g.write("CONNECT_OPTIONS={}\n".format(connections_opt))
        g.write("XOVER_FACTOR={}\n".format(crossover_factor))
        g.write("XOVER_MODIFIER={}\n".format(offset_crossover_density))
        g.write("SELECTED SEQUENCES={}\n".format(opt_console_setseq))
        g.write("SEQUENCE ORDER={}\n".format(config.AVAIL_SEQUENCES))
        g.write("PARAMETERS={}\n".format(paramstring))
        g.write("XOVER THRESHOLD={}\n".format(config.VALIDXOVERTHRESHBP))
        g.write("XOVER ADJ HELIX SPACING={}\n".format(config.VALIDXOVERSPACING_ADJ))
        g.write("XOVER SAME HELIX SPACING={}\n".format(config.VALIDXOVERSPACING_SAME))
        if os.path.isfile(os.path.join(output_dir, filename + '_seq.csv')):
            g.write("HASH={}\n".format(log.hash_file(os.path.join(output_dir, filename + '_seq.csv'))))
        g.close()
        config.save(output_dir)
    _email_success()
    job_ended = datetime.datetime.now()
    job_duration = job_ended - job_started
    log.system("Job finished in {} seconds.".format(job_duration.seconds))
    f.close()
    return 'done'
