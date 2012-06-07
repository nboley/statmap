#!/usr/bin/python

import sys
import os

# add python_lib to sys.path (at the front, since there's a trace module in 
# the standard library )
sys.path.insert(0, "../python_lib/" )

from trace import *
from wiggle import *

def usage():
    print "Usage: ./convert_trace_into_wiggle input.bin.trace [filter_threshold] > output.wig"
    sys.exit(1)

def main():
    if len(sys.argv) not in (2, 3): usage()

    # load the trace into memory
    trace = load_c_trace_from_file( sys.argv[1] )

    filter_threshold = 0
    if len(sys.argv) == 3:
        filter_threshold = float( sys.argv[2] )

    write_wiggle_from_trace_to_stream( trace, sys.stdout, filter_threshold )

if __name__ == "__main__": main()
