import shutil
import subprocess
import os
import tenkit.supernova as tk_sn
import martian

def split(args):
    mem_gb = 2048
    threads = 28
    chunk_defs = [{'__mem_gb': mem_gb, '__threads':threads, '__special': 'asmlarge' }]
    return {'chunks': chunk_defs}

def join(args, outs, chunk_defs, chunk_outs):
    shutil.move(chunk_outs[0].default, outs.default)
    summary_cs = os.path.join(outs.default, "stats", "summary_cs.csv")
    report     = os.path.join(outs.default, "stats", "summary.txt")
    shutil.copy( summary_cs, outs.summary_cs )
    shutil.copy( report, outs.report )

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
    return msg

def main(args, outs):
    print "__threads=",args.__threads
    print "__mem_gb=",args.__mem_gb
    input_dir = os.path.join(args.parent_dir,"a.base")

    cp_command = ['CP','DIR='+input_dir, 'MAX_MEM_GB='+str(args.__mem_gb), 'NUM_THREADS='+str(args.__threads)]

    if args.known_sample_id is not None:
        cp_command.append( "SAMPLE={}".format(args.known_sample_id) )
    if args.sample_id is not None:
        cp_command.append( "CS_SAMPLE_ID={}".format(args.sample_id) )
    if args.sample_desc is not None:
        cp_command.append( "CS_SAMPLE_DESC={}".format(args.sample_desc) )
    if args.addin is not None and "CP" in args.addin:
        cp_command.extend(args.addin["CP"].split())
    if args.nodebugmem is None or not args.nodebugmem:
        cp_command.append( "TRACK_SOME_MEMORY=True" )


    ## write alerts
    tk_sn.write_stage_alerts("cp", path=args.parent_dir)

    alarm_bell = tk_sn.SupernovaAlarms(base_dir=args.parent_dir)
    try:
        subprocess.check_call( cp_command )
    except subprocess.CalledProcessError as e:
        alarm_bell.post()
        ## if we actually reach this point
        ## it means that Martian::exit() was not called from C++.
        
        ## log dmesg
        log_dmesg( )

        ## detect return code
        exit_msg = process_return_code( e.returncode )
        alarm_bell.exit(exit_msg)
 
    ## post any warnings from this stage here
    alarm_bell.post()

    ## make plots
    tk_sn.try_plot_molecule_hist(args)
    tk_sn.try_plot_kmer_spectrum(args)
    
    shutil.move( args.parent_dir, outs.default )
