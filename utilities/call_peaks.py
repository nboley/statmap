#!/usr/bin/python

import sys
import os

# add python_lib to sys.path (at the front, since there's a trace module in 
# the standard library )
sys.path.insert(0, "../python_lib/" )
from peak_calling import *
from utils import *

def usage():
    print "Usage: ./call_peaks output_directory [ sample_num ]"
    sys.exit(1)

def main():
    if len(sys.argv) not in (2, 3): usage()

    smo = StatmapOutput( sys.argv[1] )

    # if we didn't enter a sample number, assume we mean all samples
    if len(sys.argv) == 2:
        call_peaks( smo.genome )
    else:
        sample_num = int( sys.argv[2] )
        sample_fname = "sample%i" % sample_num

        call_peaks_at_local_maxima( smo.genome, sample_fname )

if __name__ == "__main__":
    main()
