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
from ctypes import *

### CONSTANTS - from config.h ###
MAX_NUM_EM_ITERATIONS           = 500
MAX_PRB_CHANGE_FOR_CONVERGENCE  = 1e-2
LHD_RATIO_STOP_VALUE            = 1.50

def usage():
    print "Usage: ./sample_mapping.py output_directory ouput_fname [num_threads]"
    sys.exit(1)

def parse_arguments():
    if len(sys.argv) not in (3,4): usage()

    if len(sys.argv) == 4:  # if the number of threads to use is given
        num_threads = int( sys.argv[3] )
    else:                   # otherwise let StatmapOutput determine it
        num_threads = 0

    ofname = os.path.abspath( sys.argv[2] )

    return sys.argv[1], ofname, num_threads

def main():
    op_dir, ofname, num_threads = parse_arguments()
    
    smo = StatmapOutput( op_dir,
                         num_threads=num_threads,
                         load_mapped_reads=True,
                         load_raw_reads=True )
    
    """
    # init the starting trace 
    track_names = ["open_chromatin_density",]
    #starting_trace = load_c_trace_from_file( 
    #    "/srv/scratch/nboley/CELEGANS_ATACSeq/TEST_win1000/MLE.trace" )
    #starting_trace = init_full_trace( smo.genome, track_names )
    starting_trace = statmap_o.init_binned_trace( 
        genome, &traces, num_tracks, track_names, c_int(100) );
    statmap_o.set_trace_to_uniform(starting_trace, c_double(100));
    #statmap_o.set_trace_to_uniform( starting_trace, c_double(1.0) )
    #statmap_o.normalize_traces( starting_trace )
    statmap_o.update_mapped_reads_from_trace(
        smo.mpd_rdb,
        smo.cond_prbs_db, 
        starting_trace, 
        statmap_o.update_ATACSeq_mapped_read_prbs
    )

    statmap_o.set_bin_and_window_size(c_int(1), c_int(500))
    full_trace = init_full_trace( smo.genome, track_names )
    statmap_o.set_trace_to_uniform( full_trace, c_double(1.0) )
    update_traces_from_mapped_reads(
        smo.mpd_rdb,
        smo.cond_prbs_db, 
        full_trace,
        statmap_o.update_ATACSeq_trace_expectation_from_location
    )
    write_c_trace_to_file( full_trace, ofname )
    return
    """
    # determine the update functions
    #assert smo.config.contents.assay_type == CAGE
    update_trace_exp_from_loc = statmap_o.update_ATACSeq_trace_expectation_from_location
    update_mapped_rd_prbs = statmap_o.update_ATACSeq_mapped_read_prbs

    print "Initializing trace"
    track_names = ["open_chromatin_density",]
    bin_size = 100
    starting_trace = init_binned_trace( 
        smo.genome, track_names, 100 );
    print "Setting trace to uniform"
    statmap_o.set_trace_to_uniform( starting_trace, c_double(float(bin_size)) )
    
    # update the mapping
    print "Updating mapping"
    statmap_o.update_mapping(
        smo.mpd_rdb,
        smo.cond_prbs_db,
        starting_trace,
        c_int( MAX_NUM_EM_ITERATIONS ),
        c_float( MAX_PRB_CHANGE_FOR_CONVERGENCE ),
        c_float( LHD_RATIO_STOP_VALUE ),
        
        update_trace_exp_from_loc,
        update_mapped_rd_prbs
    )

    close_c_traces( starting_trace )

    statmap_o.set_bin_and_window_size(c_int(1), c_int(500))
    full_trace = init_full_trace( smo.genome, track_names )
    statmap_o.set_trace_to_uniform( full_trace, c_double(1.0) )
    update_traces_from_mapped_reads(
        smo.mpd_rdb,
        smo.cond_prbs_db, 
        full_trace,
        statmap_o.update_ATACSeq_trace_expectation_from_location
    )
    write_c_trace_to_file( full_trace, ofname )
    
    # clsoe the trace
    close_c_traces( starting_trace )
    
    return


if __name__ == "__main__": 
    main()
