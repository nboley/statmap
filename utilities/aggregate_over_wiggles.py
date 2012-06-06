#!/usr/bin/python

import sys
import os

# add python_lib to sys.path
sys.path.insert(0, os.path.normpath( sys.path[0] + "/../python_lib") )
from wiggle import *

THRESHOLD = 1e-6
#FLT_EPSILON = 1.1920929e-07F

def usage():
    print "Usage: ./aggregate_over_wiggles (min|max|sum) file(s).wig > output.wig"
    sys.exit(1)

def determine_aggregate_type( arg ):
    """
    Determine trace aggregate function from command line argument
    """
    if arg == "max":
        agg_fn = statmap_o.wig_lines_max
    elif arg == "min":
        agg_fn = statmap_o.wig_lines_min
    elif arg == "sum":
        agg_fn = statmap_o.wig_lines_sum
    else:
        raise Exception("Invalid aggregate type in argument : %s" % arg )

    return agg_fn

def main():
    if len(sys.argv) < 3: usage()

    # determine the aggregate type
    agg_fn = determine_aggregate_type( sys.argv[1] )

    # open the wiggles
    wigs = [ open(wig_fn) for wig_fn in sys.argv[2:] ]
    aggregate_over_wiggles( wigs, len(sys.argv)-2, sys.stdout, THRESHOLD, agg_fn )

if __name__ == "__main__": main()
