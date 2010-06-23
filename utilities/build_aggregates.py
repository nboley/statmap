#!/usr/bin/python

import os
import sys
import subprocess

valid_agg_types = ("min", "max", "sum" )

def usage():
    print "Usage: python build_aggregates.py statmap_output_directory"

# aggregate over every file in the passed directory
def aggregate_over_wigs_in_dir( directory, output_fname, aggregate_type ):
    if aggregate_type not in valid_agg_types:
        raise ValueError, "Aggregate type '%s' must be in (%s)" \
            % (aggregate_type, ", ".join(valid_agg_types))
    
    # find the wiggles to aggregate over
    fns = [ os.path.abspath(os.path.join(directory, fn)) 
            for fn in os.listdir(directory) if fn.endswith('.wig') ]

    # build the list of argumens
    args = [ aggregate_over_wiggles_exec, aggregate_type  ]
    args.extend( fns )

    # open the output file
    f = open( output_fname, "w" )
    subprocess.call( args, stdout=f )
    f.close()

    return

def aggregate_over_bs_samples(  ):
    base = "./bootstrap_samples/all_traces/"
    bs_dirs = [ os.path.join(base, dir_name) for dir_name in os.listdir( base ) ]
    for bs_dir in bs_dirs:
        # find the number of this sample
        sample_num = int( bs_dir.split("sample")[-1] )
        for agg_type in ("min", "max"):
            output_fname = "./bootstrap_samples/%s_traces/sample%i.wig" % ( agg_type, sample_num) 
            wiggles = aggregate_over_wigs_in_dir( bs_dir, output_fname, agg_type  )

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
    aggregate_over_wigs_in_dir( "./bootstrap_samples/min_traces/", "min_traces.wig", "min")
    aggregate_over_wigs_in_dir( "./bootstrap_samples/max_traces/", "max_traces.wig", "max")
