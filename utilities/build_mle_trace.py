#!/usr/bin/python

import sys
import os

# add python_lib to sys.path
sys.path.insert(0, os.path.normpath( sys.path[0] + "/../python_lib") )

from config_parsing import *
from genome import *
from mapped_read import *
from trace import *
from utils import *

### CONSTANTS - from config.h ###
MAX_NUM_EM_ITERATIONS           = 500
MAX_PRB_CHANGE_FOR_CONVERGENCE  = 1e-2
LHD_RATIO_STOP_VALUE            = 1.50

def usage():
    print "Usage: ./sample_mapping.py output_directory ouput_fname [num_threads]"
    sys.exit(1)

def main():
    if len(sys.argv) not in (3,4): usage()

    if len(sys.argv) == 3:  # if the number of threads to use is given
        num_threads = int( sys.argv[3] )
    else:                   # otherwise let StatmapOutput determine it
        num_threads = 0
    
    smo = StatmapOutput( sys.argv[1], num_threads,
                         load_mapped_reads=True,
                         load_raw_reads=True,
                         load_nc=True )

    # init the trace 
    track_names = ["fwd_strnd_read_density", "rev_strnd_read_density"]
    starting_trace = init_trace( smo.genome, track_names )

    # init cond prbs db
    cond_prbs_db = init_cond_prbs_db_from_mpd_rdb( smo.mpd_rdb );
        
    # determine the update functions
    update_trace_exp_from_loc = statmap_o.update_CAGE_trace_expectation_from_location
    update_mapped_rd_prbs = statmap_o.update_CAGE_mapped_read_prbs
    
    # initialize the starting trace
    statmap_o.set_trace_to_uniform( starting_trace, 1.0 )
    
    # update the mapping
    statmap_o.update_mapping(
        smo.mpd_rdb,
        cond_prbs_db,
        starting_trace,
        MAX_NUM_EM_ITERATIONS,
        MAX_PRB_CHANGE_FOR_CONVERGENCE,
        LHD_RATIO_STOP_VALUE,
        
        update_trace_exp_from_loc,
        update_mapped_rd_prbs
    )
    
    # clsoe the trace
    close_c_traces( starting_trace )

if __name__ == "__main__": 
    main()
