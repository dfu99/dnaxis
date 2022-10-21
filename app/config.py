import configparser


def save(outputdir):
    """Saves the settings
    Generates a text file with all configuration settings saved"""
    import os
    with open(os.path.join(outputdir, "settings"), "a") as f:
        headers = parser.items()
        for h in headers:
            f.write("\n[{}]\n".format(h[0]))
            vals = parser.items(h[0])
            for v in vals:
                f.write("{}={}\n".format(v[0].upper(), v[1]))
    f.close()


def parse_list_from_ini(li):
    _li = li.replace(" ", "")
    _li = _li.split(",")
    if all(s.isalnum() for s in _li):
        return _li
    else:
        raise TypeError


def return_by_name():
    pass


parser = configparser.ConfigParser()
parser.read('config.ini')
MAIL_PASSWORD = parser["GENERAL"]["MAIL_PASSWORD"]

# [GENERAL]
DEBUG_MSGS = bool(int(parser["GENERAL"]["DEBUG_MSGS"]))
AVAIL_SEQUENCES = parse_list_from_ini(parser["GENERAL"]["AVAIL_SEQUENCES"])
CIRCUMFERENCE_ROUNDING = bool(int(parser["GENERAL"]["CIRCUMFERENCE_ROUNDING"]))
LOG_DEBUG = bool(int(parser["GENERAL"]["LOG_DEBUG"]))
LOG_CONSOLE = bool(int(parser["GENERAL"]["LOG_CONSOLE"]))
LOG_DEV = bool(int(parser["GENERAL"]["LOG_DEV"]))
LOG_LOG = bool(int(parser["GENERAL"]["LOG_LOG"]))

# [DIRECTORY]
JOBS_DIR = parser["DIRECTORY"]["JOBS_DIR"]
DOWNLOAD_URL = parser["DIRECTORY"]["DOWNLOAD_URL"]
SEQUENCESDIR = parser["DIRECTORY"]["SEQUENCESDIR"]

# [NICKING]
LENLOW = int(parser["NICKING"]["LENLOW"])
LENUP = int(parser["NICKING"]["LENUP"])
SCAF_NICKING = parser["NICKING"]["SCAF_NICKING"]
NICKSPACING = int(parser["NICKING"]["NICKSPACING"])
PROTECT = int(parser["NICKING"]["PROTECT"])
SCPROTECT = int(parser["NICKING"]["SCPROTECT"])
LOOP_CHECK = bool(int(parser["NICKING"]["LOOP_CHECK"]))
LOOP_CHECK_RANGE = int(parser["NICKING"]["LOOP_CHECK_RANGE"])

# [CROSSOVER]
MINTPX = int(parser["CROSSOVER"]["MINTPX"])
MAXTPX = int(parser["CROSSOVER"]["MAXTPX"])
MINBPT = int(parser["CROSSOVER"]["MINBPT"])
MAXBPT = int(parser["CROSSOVER"]["MAXBPT"])
VALIDXOVERTHRESHBP = float(parser["CROSSOVER"]["VALIDXOVERTHRESHBP"])
VALIDXOVERSPACING_SAME = int(parser["CROSSOVER"]["VALIDXOVERSPACING_SAME"])
VALIDXOVERSPACING_ADJ = int(parser["CROSSOVER"]["VALIDXOVERSPACING_ADJ"])
XOVER_FACTOR = int(parser["CROSSOVER"]["XOVER_FACTOR"])
XOVERMODE = parser["CROSSOVER"]["XOVERMODE"]
GAP_SPACING_XOVER = int(parser["CROSSOVER"]["GAP_SPACING_XOVER"])
LIMSPACINGOFFSET = int(parser["CROSSOVER"]["LIMSPACINGOFFSET"])
XMDYN_MINSPACEOFF = int(parser["CROSSOVER"]["XMDYN_MINSPACEOFF"])
THRESH_ADD = float(parser["CROSSOVER"]["THRESH_ADD"])

# [EXPERIMENTAL]
ROUTING = parser["EXPERIMENTAL"]["ROUTING"]
CALIBRATIONMODE = parser["EXPERIMENTAL"]["CALIBRATIONMODE"]
DYNAMIC_SPACING = parser["EXPERIMENTAL"]["DYNAMIC_SPACING"]
TWIST_NORMALIZED = bool(int(parser["EXPERIMENTAL"]["TWIST_NORMALIZED"]))
SEEDLEN = int(parser["EXPERIMENTAL"]["SEEDLEN"])
SEED_BORROW_LIMIT = int(parser["EXPERIMENTAL"]["SEED_BORROW_LIMIT"])

# [HELICAL GEOMETRY]
AXIAL_RISE = float(parser["HELICAL GEOMETRY"]["AXIAL_RISE"])
INTERHELICAL = float(parser["HELICAL GEOMETRY"]["INTERHELICAL"])
NUCLLAGANGLE = int(parser["HELICAL GEOMETRY"]["NUCLLAGANGLE"])
BDNA_R = float(parser["HELICAL GEOMETRY"]["BDNA_R"])

# [SIMULATED ANNEALING]
HEUXOVER_RELAX_THRESH = int(parser["SIMULATED ANNEALING"]["HEUXOVER_RELAX_THRESH"])
SA_RATE = float(parser["SIMULATED ANNEALING"]["SA_RATE"])
SA_INIT_TEMP = int(parser["SIMULATED ANNEALING"]["SA_INIT_TEMP"])
SA_ACCEPT = int(parser["SIMULATED ANNEALING"]["SA_ACCEPT"])
# Number of iterations at single temperature until giving up
SA_FAIL_COUNT = int(parser["SIMULATED ANNEALING"]["SA_FAIL_COUNT"])
