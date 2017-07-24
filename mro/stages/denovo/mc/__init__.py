import os
import glob
import shutil
import subprocess
import tenkit.supernova.alerts as alerts
import martian

def split(args):
    mem_gb = 2048
    threads = 28
    chunk_defs = [{'__mem_gb': mem_gb, '__threads':threads, '__special': 'asmlarge' }]
    return {'chunks': chunk_defs}

def join(args, outs, chunk_defs, chunk_outs):
    shutil.move(chunk_outs[0].default, outs.default)

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
    input_dir = args.parent_dir

    mc_command = ['MC','DIR='+input_dir, 'MAX_MEM_GB='+str(args.__mem_gb), 'NUM_THREADS='+str(args.__threads)]

    if args.known_sample_id is not None:
        mc_command.append( "SAMPLE={}".format(args.known_sample_id) )
    if args.addin is not None and "MC" in args.addin:
        mc_command.extend(args.addin["MC"].split())


    ## write alerts
    alerts.write_stage_alerts("mc", path=args.parent_dir)

    alarm_bell = alerts.SupernovaAlarms(base_dir=args.parent_dir)
    try:
        subprocess.check_call( mc_command )
    except subprocess.CalledProcessError as e:
        alarm_bell.post()
        
        ## log dmesg
        log_dmesg( )
        
        ## log ps output
        log_ps( )
        
        ## detect return code
        exit_msg = process_return_code( e.returncode )
        alarm_bell.exit(exit_msg)
        
    ## post any warnings from this stage here
    alarm_bell.post()

    shutil.move( args.parent_dir, outs.default )
    for f in glob.glob("*.mm"):
        shutil.move(f, os.path.join(outs.default,"stats") )
