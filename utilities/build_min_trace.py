"""
Experimental rewrite of build_min_trace in Python, using python_lib ctypes interace
"""

# assume this is in utilities
import sys
import os

# add python_lib to sys.path (at the front, since there's a trace module in 
# the standard library )
sys.path.insert(0, "../python_lib/" )

from config_parsing import *
from genome import *
from mapped_read import *
from trace import *

FL_DIST_FNAME = "estimated_fl_dist.txt"

def determine_aggregate_fn( arg ):
    """
    Return pointer to the trace aggregate function specified by arg
    """
    agg_fn_ptr = c_void_p()
    if arg == "max":
        agg_fn_ptr = statmap_o.trace_max_agg
    elif arg == "min":
        agg_fn_ptr = statmap_o.trace_min_agg
    elif arg == "sum":
        agg_fn_ptr = statmap_o.trace_sum_agg
    else:
        raise Exception("Unrecognized aggregate type in argument : %s" % arg )

    return agg_fn_ptr

def determine_update_fn( assay_type ):
    """
    Return pointer to the update trace expection function for a specific assay_type
    """
    update_fn_ptr = c_void_p()

    if assay_type == CAGE:
        update_fn_ptr = statmap_o.update_CAGE_trace_expectation_from_location
    elif assay_type == CHIP_SEQ:
        update_fn_ptr = statmap_o.update_chipseq_trace_expectation_from_location
    else:
        raise Exception("Unrecognized assay type from config: %i" % assay_type)

    return update_fn_ptr

# parse arguments
USAGE = "Usage: ./build_min_trace (min|max|sum) output.bin.trace statmap_output_directory/ sample_number [num_bootstrap_samples] [--use-ip | --use-nc]"

if len(sys.argv) not in (5, 6, 7):
    print USAGE; sys.exit(1)

( agg_type, output_fname, output_dir, sample_number ) = \
    [ arg for arg in sys.argv[1:5] ]

agg_fn = determine_aggregate_fn( agg_type )

# Set num_bootstrap_samples to specified, or use default
#num_bootstrap_samples = None
#if len(sys.argv) == 6:
#    num_bootstrap_samples = int(sys.argv[5])
num_bootstrap_samples = 25 # NUM_BOOTSTRAP_SAMPLES

# Set use_nc or use_ip if option given
use_ip = False
use_nc = False
if "--use-ip" in sys.argv: use_ip = True
if "--use-nc" in sys.argv: use_nc = True
assert ( use_ip != use_nc ) or ( use_ip == use_nc == False )

# change working directory to statmap directory
try:
    os.chdir( sys.argv[3] )
except IOError as e:
    print "Could not change to specified directory : %s" % sys.argv[3]; raise

config = load_config_from_file( "config.dat" ) # CONFIG_FILENAME?

# set global variables - IMPORTANT
# threads should be from cmd line arg - def to 1 for now
statmap_o.set_num_threads(1)
# get from config
statmap_o.set_min_num_hq_bps( config.contents.min_num_hq_bps )

# determine assay type from configuration file
assay_type = config.contents.assay_type
assert assay_type == CAGE or assay_type == CHIP_SEQ
update_trace_expectation_from_location = determine_update_fn( assay_type )

# load data structures
genome = load_genome_from_disk( "genome.bin" ) # GENOME_FNAME from config.h

mpd_rdb = open_mapped_reads_db( "mapped_reads.db" ) # MAPPED_READS_DB_FNAME
mmap_mapped_reads_db( mpd_rdb )
index_mapped_reads_db( mpd_rdb )

# load fl dist if it exists
if os.path.isfile( FL_DIST_FNAME ):
    build_fl_dist_from_filename( mpd_rdb, FL_DIST_FNAME )
    build_chipseq_bs_density( mpd_rdb )

cond_prbs_db = init_cond_prbs_db_from_mpd_rdb( mpd_rdb )

# set the trace track names and determine the file names of the traces
# we're interested in

default_track_names = [ "fwd_strnd_read_density", "rev_strnd_read_density" ]
ip_track_names = [ "IP_fwd_strand", "IP_bkwd_strand" ]
nc_track_names = [ "NC_fwd_strand", "NC_bkwd_strand" ]

track_names = None
# if this is an experiment with negative control
if config.contents.NC_rdb:
    # load the set of reads depending on command line argument
    if use_ip:
        update_trace_fname = "./samples/sample{0}.ip.bin.trace".format( sample_number )
        track_names = ip_track_names
    elif use_nc:
        update_trace_fname = "./samples/sample{0}.nc.bin.trace".format( sample_number )
        track_names = nc_track_names
    else:
        raise Exception("This experiment used a negative control - you need to specify which reads to use")
else:
    update_trace_fname = "./samples/sample{0}.bin.trace".format( sample_number )
    track_names = default_track_names

# build the first trace from the specified sample; use as the trace to update
update_trace = load_c_trace_from_file( update_trace_fname )
update_cond_prbs_from_trace_and_assay_type(
    mpd_rdb, cond_prbs_db, update_trace, genome, assay_type )

print "Bootstrapping %i samples" % num_bootstrap_samples

# build aggregate, two traces at a time
curr_trace = init_trace( genome, 2, track_names )

for i in range(1, num_bootstrap_samples):
    # current_trace is reset (with zero_trace) within this function
    bootstrap_trace( curr_trace, mpd_rdb, cond_prbs_db, update_trace_expectation_from_location )
    # aggregate bootstrapped trace into update_trace
    aggregate_over_trace_pairs( update_trace, curr_trace, agg_fn )

close_c_traces( curr_trace )
write_c_trace_to_file( update_trace, output_fname )
close_c_traces( update_trace )
