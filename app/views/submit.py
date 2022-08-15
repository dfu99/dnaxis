from app import app, config
from app.routing.sequence import setscafs
from app.routing.helper import log, mymath

from flask import session, render_template
from flask_mail import Message
from app.shapeutils import custom
import os
import zipfile
from . import mail
from app.driver import driver
import traceback


def get_string_options(s, p):
    if p in s:
        return True
    else:
        return False


def compress_results():
    outputdir = session['wdir']
    try:
        # save outputs for download
        output_file = os.path.join(outputdir, "export.zip")
        if os.path.exists(output_file):
            os.remove(output_file)
        # files = os.listdir(outputdir)  # This saves all files
        files = ['prova.top', 'prova.conf', 'sequences.csv', 'modules.csv']  # Avoids sending too many under the hood files
        with zipfile.ZipFile(output_file, 'w') as zipf:
            for f in files:
                zipf.write(os.path.join(outputdir, f), f)
        zipf.close()
    except:
        raise RuntimeError


def update_params(params):
    """
    For multiple parameter relaxation
    Goes sequentially through the dict of parameters
    Iterates with first-to-last priority
    :param params:
    :return:
    """
    for key in params:
        min = params[key][0]
        max = params[key][1]
        step = params[key][2]
        precision = params[key][3]

        current_value = getattr(config, key)

        if current_value >= max:
            pass
        elif current_value < min:
            setattr(config, key, min)
        else:
            setattr(config, key, round(current_value + step, precision))
            for key2 in params:
                if key2 == key:
                    break
                else:
                    setattr(config, key2, params[key2][0])
            break
    return key


def email_error():
    output_dir = os.path.join(session['wdir'])
    jobid = output_dir[5:]
    if session['user-email']:
        msg = Message('[DNAxiS] job failed', sender='dnaxis.webmaster@gmail.com', recipients=[session['user-email']])
        msg.body = "Rendering of your DNA nanostructure encountered a problem.\n" \
                   "Please contact dfu@cs.duke.edu to submit a bug or if you need assistance." \
                   "Your job number was {}.".format(jobid)
        mail.send(msg)
    return render_template("error.html", jobno=jobid)


def email_success():
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


# call the routing module to process the input
@app.route('/process')
def upload_process():
    import json

    filename = 'test'

    f = open(os.path.join(session['wdir'], "generator_log"), "w")
    # center padding
    pad = config.AXIAL_RISE  # Pad is only for asymmetric structures
    # Create working directory
    output_dir = os.path.join(session['wdir'])
    f.close()

    # Initiate logging of output and settings
    log.new('blank', output_dir,
            console=config.LOG_CONSOLE, debug=config.LOG_DEBUG,
            developermode=config.LOG_DEV, log=config.LOG_LOG)
    print("\n")
    log.system("Output to : {}".format(output_dir))

    # Save the pathway and connections to files
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

    # build the shape class from input
    rings = session['ringdata']
    shape_class = custom.CustomInput(rings, filename)

    # Get shape parameters
    connections_opt = (3.0, 0.0)
    crossover_factor = config.XOVER_FACTOR
    offset_crossover_density = 1
    opt_console_setseq = ['m13mp18']
    opt_console_scaf_nicking = 'asymmetric_single'
    paramstring = session['paramstring']

    # Set config options if not default
    setscafs(*opt_console_setseq)
    config.SCAF_NICKING = opt_console_scaf_nicking

    # Set iteration parameters
    iter_params = {'VALIDXOVERTHRESHBP': [3.0, 2.7, -0.3, mymath.get_precision(0.1)]
                   }

    # Initialize
    for iter_key in iter_params:
        setattr(config, iter_key, iter_params[iter_key][0])

    iterate_all = False  # Keep going even if successful?
    active_key = list(iter_params.keys())[0]

    log.system("\nSet {} = {}".format(active_key, getattr(config, active_key)))
    try:
        driver(filename,
               output_dir,
               shape_class,
               crossover_factor,
               connections_opt,
               offset_crossover_density,
               auto_scaf_options=get_string_options(paramstring, '-aso'),
               shape_only=get_string_options(paramstring, '-shape'),
               skip_routing=get_string_options(paramstring, '-skroute'),
               skip_nicks=get_string_options(paramstring, '-sknicks'),
               skip_sequence=get_string_options(paramstring, '-skseq'),
               force_shuffle_pathway=get_string_options(paramstring, '-fsp'),
               force_rand_seq=get_string_options(paramstring, '-frs'),
               force_reseeding=get_string_options(paramstring, '-freseed'),
               optimize_xovers=get_string_options(paramstring, '-ox'),
               use_extensions=get_string_options(paramstring, '-ext'),
               add_uvxlinking=get_string_options(paramstring, '-uvxl'),
               force_xover_density=get_string_options(paramstring, '-fxd'),
               replace_long_bonds=get_string_options(paramstring, '-longbonds'),
               enable_validate=get_string_options(paramstring, '-valid'),
               enable_stats=get_string_options(paramstring, '-stats'),
               enable_debugger=get_string_options(paramstring, '-debug'),
               force_debug_procedures=get_string_options(paramstring, '-debug'),
               force_clear_seam=get_string_options(paramstring, '-fcs'),
               force_clean_merge=get_string_options(paramstring, '-fcm'),
               use_old_routing=get_string_options(paramstring, '-oldrouting'),
               twist_normalized=get_string_options(paramstring, '-tn'),
               save_steps=get_string_options(paramstring, '-savesteps'))

        # Save shape specific settings to file
        g = open(os.path.join(output_dir, "settings"), "w")
        g.write("[SHAPE]\n")
        g.write("NAME={}\n".format(filename))
        for param in shape_class.inputs:
            g.write("{}={}\n".format(param, shape_class.inputs[param]))
        g.write("CONNECT_OPTIONS={}\n".format(connections_opt))
        g.write("XOVER_FACTOR={}\n".format(crossover_factor))
        g.write("XOVER_MODIFIER={}\n".format(offset_crossover_density))
        g.write("SELECTED SEQUENCES={}\n".format(opt_console_setseq))
        g.write("PARAMETERS={}\n".format(paramstring))
        if os.path.isfile(os.path.join(output_dir, filename + '_seq.csv')):
            g.write("HASH={}\n".format(log.hash_file(os.path.join(output_dir, filename + '_seq.csv'))))
        g.close()
        config.save(output_dir)
        # Compress results into ZIP file
        try:
            compress_results()
        except RuntimeError:
            log.system("Could not collect all the files for export.zip.")
            f.close()
            raise RuntimeError
    except Exception as e:
        email_error()
        with open(os.path.join(output_dir, "error.txt"), 'a') as f:
            f.write(str(e))
            f.write(traceback.format_exc())
        log.system("Encountered an unknown error. If user supplied an email, they were notified.")
        f.close()
    email_success()

    f.close()
