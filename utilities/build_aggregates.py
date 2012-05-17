#!/usr/bin/python
"""
Build aggregate traces from the traces in ./bootstrap_samples/all_traces
using aggregate_over_traces

For performance, you can specify a single aggregation type to use as an
optional argument. The default is to generate max and min aggregates.
"""

import os
import sys
import subprocess

# set in main
aggregate_exec = None
assay_type = None
agg_types = None

valid_agg_types = ("min", "max", "sum" )
assay_types = { # mapping integer enum values of assay_type_t to assays
    1 : 'cage',
    2 : 'chipseq',
}

def usage():
    print "Usage: python aggregate_over_all_traces.py statmap_output_directory [%s]" % (
            '|'.join( valid_agg_types )
        )

# aggregate over every file in the passed directory
def aggregate_over_traces_in_dir( directory, output_fname, aggregate_type, suffix=".bin.trace" ):
    """Aggregate over all of the traces in a diretory.
    
    Suffix is a bit of a hack. I sometimes need to aggregate over *just* the nc
    ( negative control ) or ip ( immino-precipitate ) traces - this provides
    a ( very clunk ) mechanism. It defaults to empty.
    """
    if aggregate_type not in valid_agg_types:
        raise ValueError, "Aggregate type '%s' must be in (%s)" \
            % (aggregate_type, ", ".join(valid_agg_types))
    
    # find the wiggles to aggregate over
    fns = [ os.path.abspath(os.path.join(directory, fn)) 
            for fn in os.listdir(directory) if fn.endswith( suffix ) ]

    # build the list of arguments
    args = [ aggregate_exec, aggregate_type, output_fname  ]
    args.extend( fns )
    
    try: 
        rv = subprocess.check_call( args, stdout=sys.stdout )
    except subprocess.CalledProcessError:
        print args
        print " ".join( args )
        raise
    
    return

def aggregate_over_bs_samples(  ):
    """
    Loop over the bootstrap traces in bootstrap_samples/all_traces to generate
    bootstrap_samples/max_traces and bootstrap_samples/min_traces
    """
    base = "./bootstrap_samples/all_traces/"
    bs_dirs = [ os.path.join(base, os.path.split(dir_name)[-1]) for dir_name in os.listdir( base ) ]
    for bs_dir in bs_dirs:
        # find the number of this sample
        sample_num = int( bs_dir.split("sample")[-1].split(".")[0] )
        # aggregate differently depending on assay type
        if assay_type == "chipseq":
            for track_type in ("nc", "ip" ):
                for agg_type in agg_types:
                    output_fname = "./bootstrap_samples/%s_traces/sample%i.%s.bin.trace" % (
                            agg_type, sample_num, track_type ) 
                    wiggles = aggregate_over_traces_in_dir(
                            bs_dir, output_fname, agg_type, ".%s.bin.trace" % track_type  )
        elif assay_type == "cage":
            for agg_type in agg_types:
                output_fname = "./bootstrap_samples/%s_traces/sample%i.bin.trace" % (
                        agg_type, sample_num )
                wiggles = aggregate_over_traces_in_dir(
                        bs_dir, output_fname, agg_type, ".bin.trace" )

if __name__ == "__main__":
    if len( sys.argv ) not in ( 2, 3 ):
        usage()
        sys.exit(-1)
    
    ## parse the arguments, and setup the environment
    # the directory of the statmap output
    output_dir = sys.argv[1]    
    # the path to this script. The C aggregate utility will be in the same directory.
    aggregate_exec = os.path.join( \
        os.path.abspath( os.path.split(sys.argv[0])[0] ), "aggregate_over_traces" \
    )
    os.chdir(output_dir)

    # determine assay type from config.dat (integer value of assay_type_t enum)
    with open("config.dat") as conf_fp:
        for line in conf_fp:
            if line.split()[0] == "assay_type:":
                assay_type = assay_types[ int(line.split()[1]) ]

    # set aggregate type(s)
    if len( sys.argv ) == 3: # set by cmd line argument
        agg_type = sys.argv[2]

        # make sure user input is a valid aggregate type
        if agg_type not in valid_agg_types:
            usage()
            sys.exit(-1)

        agg_types = ( agg_type, )
    else: # default
        agg_types = ( "min", "max" )

    # aggregate over the bootstrap samples
    print "Aggregating over bootstrap samples..."
    aggregate_over_bs_samples( )

    # aggregate over the aggregated bootstrap samples
    print "Aggregating over the aggregates..."
    if assay_type == "chipseq":
        for agg_type in agg_types:
            for track_type in ("nc", "ip"):
                aggregate_over_traces_in_dir(
                    "./bootstrap_samples/%s_traces/" % agg_type,
                    "%s_trace.%s.bin.trace" % ( agg_type, track_type ),
                    agg_type,
                    ".%s.bin.trace" % track_type
                )
    elif assay_type == "cage":
        for agg_type in agg_types:
            aggregate_over_traces_in_dir(
                "./bootstrap_samples/%s_traces/" % agg_type,
                "%s_trace.bin.trace" % agg_type,
                agg_type,
                ".bin.trace"
            )
