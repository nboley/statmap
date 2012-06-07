#!/usr/bin/python

import sys
import os

# add python_lib to sys.path
sys.path.insert(0, os.path.normpath( sys.path[0] + "/../python_lib") )
from trace import *
from utils import *

def usage():
    print "Usage: ./build_marginal_mappings_trace_from_mpd_read_db output_directory [filter_threshold] > output.wig"
    sys.exit(1)

def parse_args():
    if len(sys.argv) not in (2, 3): usage()

    try:
        filter_threshold = float( sys.argv[2] )
    except:
        filter_threshold = 0.0

    smo = StatmapOutput( sys.argv[1],
                         load_mapped_reads=True )

    return (smo, filter_threshold)

def main():
    smo, filter_threshold = parse_args()

    reset_all_read_cond_probs( smo.mpd_rdb, smo.cond_prbs_db )

    # initialize a trace
    track_names = [ "fwd_strand", "bkwd_strand" ]
    trace = init_trace( smo.genome, 2, track_names )

    # update the trace from the naive kernel
    # we hack this a bit by noting that the naive kernel is just the
    # trace kernel (stranded, which we nearly always want )
    update_traces_from_mapped_reads(
        smo.mpd_rdb, smo.cond_prbs_db, trace, smo.trace_update_fn
    )

    write_wiggle_from_trace_to_stream( trace, sys.stdout, filter_threshold )

if __name__ == "__main__": main()
