import os
import glob
import shutil
import subprocess
import tenkit.supernova.alerts as alerts
import martian

def split(args):
    mem_gb = 2048
    threads = 28
    chunk_defs = [{'__mem_gb': mem_gb, '__threads':threads, '__special': 'asmlarge'}]
    return {'chunks': chunk_defs}

def join(args, outs, chunk_defs, chunk_outs):
    shutil.move(chunk_outs[0].default, outs.default)

def check_exclude(path, ext):
    head=path[:-len(ext)]
    tail=path[-len(ext):]
    if tail != ext:
        raise Exception("file has incorrect extension: " + path)
    return head

def probe_reads(path):
    out=subprocess.check_output(['FastFastbCount','QUIET=True','FASTB='+path]).strip()
    return long(out)

def probe_read_len(path):
    # below is harcoded to 150 for internal use only 
    # CS pipeline uses target_reads path only
    return 150

def log_dmesg( ):
    num_lines = 100
    try:
        dmesg_out = subprocess.check_output( ["dmesg"] )
        martian.log_info( "----dmesg output----" )
        dmesg_lines = dmesg_out.split("\n")
        if len(dmesg_lines) <= num_lines:
            martian.log_info( dmesg_out )
        else:
            for line in dmesg_lines[-num_lines:]:
                martian.log_info( line )
        martian.log_info( "--------------------" )
    except Exception as e:
        martian.log_info( "Unable to run dmesg." )
        martian.log_info( str(e) )

def log_ps( ):
    MAX_LINES = 6
    ps_cmd = ["ps", "--sort=-rss", "-eo", "pid,pmem,rss,comm,uid"]
    ps_log = "Running command " + " ".join(ps_cmd) + "\n"
    try:
        ps_out = subprocess.check_output( ps_cmd )
        for i, line in enumerate( ps_out.split("\n") ):
            if i == MAX_LINES:
                break
            ps_log += "%8s %5s %12s %5.5s %s\n" % tuple(line.split())
    except:
        ps_log += "Running ps command failed\n"
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

def main(args, outs):
    print "__threads=",args.__threads
    print "__mem_gb=",args.__mem_gb
    h1 = check_exclude(args.reads, '.fastb')
    h2 = check_exclude(args.quals, '.qualp')
    h3 = check_exclude(args.bci, '.bci')
    if h1 != h2 or h2 != h3:
        raise Exception( "something wrong with filenames passed in" )

    select_frac=1.0
    if args.downsample is not None:
        if args.downsample.get("target_reads", None) is not None:
            target_nreads = args.downsample["target_reads"]
            actual_nreads = probe_reads( args.reads )

            print "target_nreads={}, actual_nreads={}".format( \
                target_nreads, actual_nreads )

            if target_nreads > actual_nreads:
                select_frac=1.
            else:
                select_frac=1.*target_nreads/actual_nreads

        elif args.downsample.get("gigabases", None) is not None:
            target_gigabases = args.downsample["gigabases"]
            actual_reads = probe_reads(args.reads)
            actual_maxlen = probe_read_len(args.reads)
            actual_gigabases = actual_reads*actual_maxlen/1.e9

            print "target_gigabases={}, actual_gigabases={}, actual_reads={}, actual_maxlen={}".format( \
                target_gigabases, actual_gigabases, actual_reads, actual_maxlen )

            if target_gigabases > actual_gigabases:
                select_frac=1.
            else:
                select_frac=1.*target_gigabases/actual_gigabases

        else:
            martian.log_warn("unknown downsample mode was ignored")

    select_frac=1.0

    df_command = ['DF', 'LR_SELECT_FRAC={:.8f}'.format(select_frac), 'LR='+args.reads,
            'OUT_DIR='+outs.default, "MAX_MEM_GB="+str(args.__mem_gb), "NUM_THREADS="+str(args.__threads) ]

    if args.mspedges is not None:
        df_command.append( "MSPEDGES={}".format(args.mspedges) )

    if args.pipeline_id is not None:
        df_command.append( "PIPELINE={}".format(args.pipeline_id) )

    if args.known_sample_id is not None:
        df_command.append( "SAMPLE={}".format(args.known_sample_id) )

    if args.addin is not None and "DF" in args.addin:
        df_command.extend(args.addin["DF"].split())

    if args.nodebugmem is None or not args.nodebugmem:
        df_command.append( "TRACK_SOME_MEMORY=True" )

    print " ".join(df_command)
    
    ## write alerts to txt file
    alerts.write_stage_alerts("df", path=outs.default)

    alarm_bell = alerts.SupernovaAlarms(base_dir=outs.default)
    try:
        subprocess.check_call(df_command)
    except subprocess.CalledProcessError as e:
        ## there was an error of some kind
        ## if any Martian::exit was issued in the C++ code
        ## then we call martian.exit from here
        alarm_bell.post()

        ## if we actually reach this point
        ## it means that Martian::exit() was not called from C++.
        
        ## log dmesg
        log_dmesg( )
        
        ## log ps output
        log_ps( )

        ## detect return code
        exit_msg = process_return_code( e.returncode )
        alarm_bell.exit(exit_msg)
    
    ## if DF completed successfully,
    ## post any warnings from this stage here
    alarm_bell.post()

    for f in glob.glob("*.mm"):
        shutil.move(f, os.path.join(outs.default,"stats") )
