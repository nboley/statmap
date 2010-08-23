#!/usr/bin/python

import os
import sys
import subprocess

valid_agg_types = ("min", "max", "sum" )

def usage():
    print "Usage: python build_aggregates.py statmap_output_directory"

# aggregate over every file in the passed directory
def aggregate_over_wigs_in_dir( directory, output_fname, aggregate_type, suffix=".wig" ):
    """Aggregate over all of the wiggles in a diretory.
    
    Suffix is a bit of a hack. I sometimes need to aggregate over *just* the nc
    ( negative control ) or ip ( immino-precipitate wiggles )wiggles - this provides
    a ( very clunk ) mechanism. It defaults to empty.  
    """
    if aggregate_type not in valid_agg_types:
        raise ValueError, "Aggregate type '%s' must be in (%s)" \
            % (aggregate_type, ", ".join(valid_agg_types))
    
    # find the wiggles to aggregate over
    fns = [ os.path.abspath(os.path.join(directory, fn)) 
            for fn in os.listdir(directory) if fn.endswith( suffix ) ]

    # build the list of argumens
    args = [ aggregate_over_wiggles_exec, aggregate_type  ]
    args.extend( fns )

    # open the output file
    f = open( output_fname, "w" )
    try: 
        rv = subprocess.check_call( args, stdout=f )
    except subprocess.CalledProcessError:
        print args
        print " ".join( args )
        raise
    f.close()

    return

def aggregate_over_bs_samples(  ):
    base = "./bootstrap_samples/all_traces/"
    bs_dirs = [ os.path.join(base, os.path.split(dir_name)[-1]) for dir_name in os.listdir( base ) ]
    for bs_dir in bs_dirs:
        # find the number of this sample
        sample_num = int( bs_dir.split("sample")[-1].split(".")[0] )
        for track_type in ("nc", "ip" ):
            for agg_type in ("min", "max"):
                output_fname = "./bootstrap_samples/%s_traces/sample%i.%s.wig" % ( agg_type, sample_num, track_type) 
                wiggles = aggregate_over_wigs_in_dir( bs_dir, output_fname, agg_type, ".%s.wig" % track_type  )

if __name__ == "__main__":
    if len( sys.argv ) != 2:
        usage()
        sys.exit(-1)
    
    ## parse the arguments, and setup the environment
    # the directory of the statmap output
    output_dir = sys.argv[1]    
    # the path to this script. The C aggregate utility will be in the same directory.
    global aggregate_over_wiggles_exec
    aggregate_over_wiggles_exec = os.path.join( \
        os.path.abspath( os.path.split(sys.argv[0])[0] ), "aggregate_over_wiggles" \
    )
    os.chdir(output_dir)
    
    # aggregate over the bootstrap samples
    aggregate_over_bs_samples( )
    # aggregate over the aggregated bootstrap samples
    aggregate_over_wigs_in_dir( "./bootstrap_samples/min_traces/", "min_trace.nc.wig", "min", ".nc.wig")
    aggregate_over_wigs_in_dir( "./bootstrap_samples/max_traces/", "max_trace.nc.wig", "max", ".nc.wig")
    aggregate_over_wigs_in_dir( "./bootstrap_samples/min_traces/", "min_trace.ip.wig", "min", ".ip.wig")
    aggregate_over_wigs_in_dir( "./bootstrap_samples/max_traces/", "max_trace.ip.wig", "max", ".ip.wig")
