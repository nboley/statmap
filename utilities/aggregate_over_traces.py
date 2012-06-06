#!/usr/bin/python

import sys
import os

# add python_lib to sys.path
sys.path.insert(0, os.path.normpath( sys.path[0] + "/../python_lib") )
from trace import *

def usage():
    print "Usage: ./aggregate_over_traces.py (min|max|sum) output.bin.trace file(s).bin.trace"
    sys.exit(1)

def determine_aggregate_type( arg ):
    """
    Determine trace aggregate function from command line argument
    """
    if arg == "max":
        agg_fn_ptr = statmap_o.trace_agg_max
    elif arg == "min":
        agg_fn_ptr = statmap_o.trace_agg_min
    elif arg == "sum":
        agg_fn_ptr = statmap_o.trace_agg_sum
    else:
        raise Exception("Invalid aggregate type: %s" % arg )

    return agg_fn_ptr

def main():
    if len(sys.argv) <= 3: usage()

    agg_fn = determine_aggregate_type( sys.argv[1] )

    # load the first trace as the trace we will update
    update_trace = load_c_trace_from_file( sys.argv[3] )

    # loop over the .bin.trace input files
    for trace_fname in sys.argv[4:]:
        curr_trace = load_c_trace_from_file( trace_fname )
        aggregate_over_trace_pairs( update_trace, curr_trace, agg_fn )
        close_c_traces( curr_trace )

    # write final update trace to file
    write_c_trace_to_file( update_trace, sys.argv[2] )
    close_c_traces( update_trace )

if __name__ == "__main__": main()
