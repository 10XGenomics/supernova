#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
# Subprocess support code for Supernova

import martian
import alerts
import textwrap
import subprocess
import shlex
import shutil
import types
import exceptions
import glob
import sys


class CxxStage:
    """aggregate some common python stage code behavior for running subprocesses
    arguments:
    name is the lowercase stage name as would be used in the addin dictionary
        (e.g. "df", "cp", "tr", etc.)
    path is the input path -- outs.default for the first in the chain,
                              args.parent_dir otherwise
    dstpath is the output path -- generally outs.default
    process is a list of things to do in the try section:
        - items that are functions are called
        - items that are strings are passed to shlex.split() and subprocess.check_call()
        - items that are lists are assumed to be args for subprocess.check_call()
    post_succeed and post_failure are both functions for after-work in the case
        that the subprocess succeeds or fails
    """
    def __init__(self, name, path, dstpath, process, \
                    post_succeed = lambda: None, post_failure = lambda: None):
        self.name = name
        self.path = path
        self.dstpath = dstpath
        self.process = process
        self.post_succeed = post_succeed
        self.post_failure = post_failure
        self.lastpid = None

        self.rollup()

    def run_subcommand(self, args):
        try:
            proc = subprocess.Popen( args, 0 )
        except Exception as e:
            raise RuntimeError( "Popen failed with args: " + str(args) +  \
                    ", Exception: "+str(e) )
        self.lastpid = proc.pid
        proc.wait()
        if proc.returncode != 0:
            raise subprocess.CalledProcessError( returncode = proc.returncode, \
                    cmd = " ".join(args), output = None)

    def rollup(self):
        ## write alerts to txt file
        alarms_and_alerts = self.path + "/alarms"
        alerts.write_stage_alerts(self.name, path=alarms_and_alerts)
        alarm_bell = alerts.SupernovaAlarms(base_dir=alarms_and_alerts)
        try:
            for p in self.process:
                if type(p) is types.ListType:
                    self.run_subcommand( p )
                elif type(p) is types.StringType:
                    self.run_subcommand(shlex.split(p))
                elif type(p) is types.FunctionType:
                    p()
                else:
                    raise exceptions.RuntimeError("unknown type process passed to rollup")

        except subprocess.CalledProcessError as e:
            ## there was an error of some kind
            ## if any Martian::exit was issued in the C++ code
            ## then we call martian.exit from here
            alarm_bell.post()

            self.post_failure()

            ## if we actually reach this point
            ## it means that Martian::exit() was not called from C++.

            ## log dmesg
            log_dmesg( self.lastpid )

            ## log ps output
            log_ps( )

            ## detect return code
            exit_msg = process_return_code( e.returncode )
            alarm_bell.exit(exit_msg)

        ## if completed successfully,
        ## post any warnings from this stage here
        alarm_bell.post()

        self.post_succeed()

        if self.dstpath is not None and self.dstpath != self.path:
            shutil.move( self.path, self.dstpath )
        for f in glob.glob("*.mm"):
            shutil.move(f, os.path.join(path,"stats") )


def log_dmesg( pid = None ):
    num_lines = 100
    try:
        output=textwrap.dedent("""
               ---- START of dmesg output ----
               This is output of the dmesg command captured by Supernova(tm).  The primary
               intent is to find notices from the system-wide Out-of-Memory killer that
               would indicate that the Supernova process was terminated by the system
               for excessive memory use.

               Note that 'dmesg' has no timestamp associated with it, so messages seen
               may be quite old.

               """)
        if pid is not None:
            output+="Please look for process id " + str(pid) + "\n\n"

        output+="\n".join( map( lambda s: "\t"+s,
                    subprocess.check_output( ["dmesg"] ).split("\n")[-num_lines:] ) )
        output+="---- END of dmesg output ----\n"
        martian.log_info(output)
    except Exception as e:
        martian.log_info( "Unable to run dmesg: " + str(e) )

def log_ps( pid = None ):
    MAX_LINES = 6
    ps_cmd = ["ps", "--sort=-rss", "-eo", "pid,pmem,rss,uid,cmd"]

    ps_log=textwrap.dedent("""
            ---- START of ps output ----
            Output of the 'ps' command is intended to diagnose whether there are
            other memory intensive processes running on the same machine.  Process
            names are limited to the first five characters, and Supernova processes
            show up as "exe."  Note that the Supernova process that failed here will
            *not* be in the list, as it has already exited.  The RSS column describes
            the physical memory used by the process in kB.

            """)

    if pid is not None:
       ps_log+="Note that he process id for this Supernova process was " + str(pid) + "\n\n"

    ps_log += "Running command " + " ".join(ps_cmd) + "\n"
    try:
        ps_out = subprocess.check_output( ps_cmd )
        for line in ps_out.split("\n")[-MAX_LINES:]:
            s=line.split()
            if len(s) == 5:
                ps_log += "\t%8s %5s %12s %s %5.5s\n" % tuple(s)
            else:
                ps_log += " ".join(s)
    except Exception as e:
        ps_log += "Running ps command failed: " + str(e) + " on line " + line + "\n"
    ps_log += "---- END of ps output ----\n"
    martian.log_info( ps_log )

def process_return_code( returncode ):
    msg = None
    if returncode < 0:
        SIG_DICT = { -9: "KILL signal",
                     -1: "HUP signal",
                     -2: "INT signal",
                     -15: "TERM signal" }
        sig_name = SIG_DICT.get( returncode, "signal" )
        msg = "A Supernova process was terminated with a %s "\
        "(code: %d). This may have been sent by you, your IT admin, "\
        "or automatically by the system itself "\
        "(e.g. the out-of-memory killer)." % (sig_name,-returncode)
    elif returncode == 99:
        msg = "Supernova terminated because of insufficient memory. "\
        "The stage _stdout file may contain additional useful information, "\
        "e.g., whether any competing processes were running."
    return msg

def standard_split(special=None, maxcores=None):
    mem_gb = 2048
    if maxcores is None: maxcores=128   # should get knocked down by Martian
    chunk_defs = [{'__mem_gb': mem_gb, '__threads':maxcores}]
    if special is not None: chunk_defs[0]['__special'] = special
    return {'chunks': chunk_defs}

def append_if( lst, value, fmt ):
    if value is not None and value != "":
        lst.append( fmt.format(value) )

def standard_args( cmdlist, args, addin_name, tsm=True ):
    cmdlist.append( "MAX_MEM_GB={}".format( args.__mem_gb ) )
    cmdlist.append( "NUM_THREADS={}".format( args.__threads ) )

    if args.addin is not None and addin_name in args.addin:
        cmdlist.extend( args.addin[addin_name].split() )

    if tsm:
        if args.nodebugmem is None or not args.nodebugmem:
            cmdlist.append( "TRACK_SOME_MEMORY=True" )
