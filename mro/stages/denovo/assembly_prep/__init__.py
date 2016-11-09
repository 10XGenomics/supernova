import subprocess
import os

def check_file_type( filename, extension ):
    if filename[ -len(extension): ] != extension:
        raise Exception("file " + filename +
                " doesn't have expected extension " + extension )
    return filename[ :-len(extension) ]


def split(args):
    chunks=[ {'in_reads_file' : file} for file in args.in_reads ]
    return { 'chunks': chunks }


def join(args, outs, chunk_defs, chunk_outs):
    out_dir = os.path.dirname( outs.out_reads )
    out_head = os.path.join(out_dir, "out_reads")

    with open( 'merge_heads.arg', 'w' ) as f:
        # okay, so our code takes a list of "heads"
        # we'll sanity check the chunk_outs and see if it matches
        for chunk_out in chunk_outs:
            h1=check_file_type( chunk_out.out_reads, '.fastb' )
            h2=check_file_type( chunk_out.out_quals, '.qualp' )
            if ( h1 != h2 ) :
                raise Exception( "join() has unaligned outs :(" )
            f.write( h1 + "\n" )

    subprocess.check_call([ 'ParseBarcodedFastqs',
                           'OUT_HEAD='+out_head,
                           'MERGE_HEADS=@merge_heads.arg' ])

    os.unlink("merge_heads.arg")

    outs.out_reads = out_head + ".fastb"
    outs.out_quals = out_head + ".qualp"
    outs.out_bci   = out_head + ".bci"
#    os.rename( out_head + ".fastb", outs.out_reads )
#    os.rename( out_head + ".qualp", outs.out_quals )


def main(args, outs):
    check_file_type(outs.out_reads, ".fastb")
    check_file_type(outs.out_quals, ".qualp")
    check_file_type(outs.out_bci, ".bci")

    out_dir = os.path.dirname( outs.out_reads )
    out_head = os.path.join(out_dir, "out_reads")

    subprocess.check_call([ 'ParseBarcodedFastqs',
                           'OUT_HEAD='+out_head,
                           'FASTQS='+args.in_reads_file ] )

    outs.out_reads = out_head + ".fastb"
    outs.out_quals = out_head + ".qualp"
    outs.out_bci   = out_head + ".qualp"
#    os.rename( out_head + ".fastb", outs.out_reads )
#    os.rename( out_head + ".qualp", outs.out_quals )
