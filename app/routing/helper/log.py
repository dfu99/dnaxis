#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 23:38:02 2018

@author: dfu

filename: log

consolidated logging function for printing messages to file

changelog:
    9/12/2018:
        added debug label to printout
        added option to add headers
"""
import datetime
import os

global __DEBUG__
__DEBUG__ = True

global __CONSOLE__
__CONSOLE__ = True

global __LOG__
__LOG__ = True

global __DEV__
__DEV__ = True

global __WDIR__
__WDIR__ = ""

global __SDIR__
__SDIR__ = ""


def new(fout, wdir, debug=True, console=True, log=True, developermode=False):
    global __DEBUG__
    __DEBUG__ = debug
    global __CONSOLE__
    __CONSOLE__ = console
    global __LOG__
    __LOG__ = log
    global __DEV__
    __DEV__ = developermode
    global filename
    filename = fout

    log_fileloc = os.path.join(wdir, 'log.txt')
    os.makedirs(os.path.dirname(log_fileloc), exist_ok=True)

    sys_fileloc = os.path.join(wdir, 'last_console_output.txt')
    os.makedirs(os.path.dirname(sys_fileloc), exist_ok=True)

    global __WDIR__
    __WDIR__ = log_fileloc

    global __SDIR__
    __SDIR__ = sys_fileloc

    with open(__WDIR__, 'w') as f:
        f.write(str(datetime.datetime.now()) + "\n")
        f.write("Output options: " +
                "stdout({}) ".format(__CONSOLE__) +
                "log ({}) ".format(__LOG__) +
                "debug({}) ".format(__DEBUG__) +
                "developer({})\n".format(__DEV__))
    f.close()

    with open(__SDIR__, 'w') as f:
        f.write(str(datetime.datetime.now()) + "\n")
    f.close()


def out(name, *msg, label="INFO", headerlevel=0):
    #    if name == "routing.debug" and not __DEBUG__:
    #        return
    #    elif name == "routing.debug" and __DEBUG__:
    #        label = "DEBUG"

    msg = [str(m) for m in msg]
    outmsg = " ".join(msg)
    if headerlevel != 0:
        outmsg = padheader(outmsg, headerlevel)
    outstr = "[{}]:{}:{}\n".format(label, name, outmsg)
    if __CONSOLE__:
        print(outstr)
    if __LOG__:
        with open(__WDIR__, 'a') as f:
            f.write(outstr)
        f.close()


def debug(name, *msg, label="DEBUG", headerlevel=0):
    if __DEBUG__:
        msg = [str(m) for m in msg]
        outmsg = " ".join(msg)
        if headerlevel != 0:
            outmsg = padheader(outmsg, headerlevel)
        outstr = "[{}]:{}:{}\n".format(label, name, outmsg)
        if __CONSOLE__:
            print(outstr)
        if __LOG__:
            with open(__WDIR__, 'a') as f:
                f.write(outstr)
            f.close()


def dev(name, *msg, label="DEVELOPER", headerlevel=0):
    if __DEV__:
        msg = [str(m) for m in msg]
        outmsg = " ".join(msg)
        if headerlevel != 0:
            outmsg = padheader(outmsg, headerlevel)
        outstr = "[{}]:{}:{}\n".format(label, name, outmsg)
        if __CONSOLE__:
            print(outstr)
        if __LOG__:
            with open(__WDIR__, 'a') as f:
                f.write(outstr)
            f.close()


def makeheader(outstr):
    padout = "{: ^100}".format(outstr)
    return padout


def padheader(outstr, headerlevel):
    if headerlevel == 1:
        padout = "{:=^80}".format(outstr)
    elif headerlevel == 2:
        padout = "{:*^70}".format(outstr)
    elif headerlevel == 3:
        padout = "{:%^60}".format(outstr)
    elif headerlevel == 4:
        padout = "{:+^50}".format(outstr)
    padout = makeheader(padout)
    return padout


def log_error(e):
    """
    Log errors that interrupt runtime
    :param e: Exception string
    :return:
    """
    print("[ERROR]:", e)
    out(__name__, "[ERROR]:", e)
    raise RuntimeError(e)


def log_warning(e):
    """
    Log warnings that don't interrupt runtime
    :param e: Exception string
    :return:
    """
    if __DEBUG__:
        print("[WARNING]:", e)
    out(__name__, "[WARNING]:", e)


def log_debug(e):
    """
    Print temporarily used debug messages
    :param e:
    :return:
    """
    if __DEBUG__:
        print("[DEBUG]:", e)
    out(__name__, "[DEBUG]:", e)


def system(m):
    """
    Print system level messages
    :param m:
    :return:
    """
    print(m)
    with open(__SDIR__, 'a') as f:
        f.write(m)
        f.write("\n")
    f.close()
    pass


def hash_file(filename):
    import hashlib
    hasher = hashlib.md5()
    if not os.path.exists(filename):
        with open(filename, 'w') as afile:
            afile.write("")
        afile.close()
    with open(filename, 'rb') as afile:
        buf = afile.read()
        hasher.update(buf)
    return hasher.hexdigest()
