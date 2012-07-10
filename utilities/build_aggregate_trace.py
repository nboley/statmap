#!/usr/bin/python

import sys
import os

# add python_lib to sys.path
sys.path.insert(0, os.path.normpath( sys.path[0] + "/../python_lib") )
from trace import *
from utils import *

def get_aggregate_fn( arg ):
    if arg == "max":
        agg_fn_ptr = statmap_o.trace_agg_max
    elif arg == "min":
        agg_fn_ptr = statmap_o.trace_agg_min
    elif arg == "sum":
        agg_fn_ptr = statmap_o.trace_agg_sum
    else:
        raise Exception("Unrecognized aggregate type in argument : %s" % arg )

    return agg_fn_ptr

def parse_args():
    import argparse; parser = argparse.ArgumentParser()

    # Required arguments
    parser.add_argument("aggregate_type",
        help="Aggregate function to use",
        choices=[ "min", "max", "sum"] )
    parser.add_argument("output_directory",
        help="Path to Statmap output directory to process")
    parser.add_argument("output_filename",
        help="Output filename for binary trace")
    parser.add_argument("sample_number",
        help="Sample number to bootstrap from",
        type=int)

    # Optional arguments
    parser.add_argument("--num-bootstrap-samples", "-nbs",
        help="Number of bootstrap samples to take. Default: %(default)s",
        type=int,
        default=25 )
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--use-ip", "-ip",
help="(For Chipseq only) Use the immunoprecipitate mapped reads (default). ",
        action="store_true")
    group.add_argument("--use-nc", "-nc",
        help="(For Chipseq only) Use the negative control mapped reads. ",
        action="store_true")

    return parser.parse_args()

def set_track_names( smo, args ):
    # define possible track names
    default_track_names = [ "fwd_strnd_read_density", "rev_strnd_read_density" ]
    ip_track_names = [ "IP_fwd_strand", "IP_bkwd_strand" ]
    nc_track_names = [ "NC_fwd_strand", "NC_bkwd_strand" ]

    # if this is a Chip-Seq experiment (has a negative control)
    if smo.config.contents.assay_type == CHIP_SEQ:
        if args.use_ip:
            update_trace_fname = \
                "./samples/sample{0}.ip.bin.trace".format( args.sample_number )
            track_names = ip_track_names
        elif args.use_nc:
            update_trace_fname = \
                "./samples/sample{0}.nc.bin.trace".format( args.sample_number )
            track_names = nc_track_names
        else:
            raise Exception(
"This experiment used a negative control - you need to specify which reads to "
"use")
    else:
        update_trace_fname = \
            "./samples/sample{0}.bin.trace".format( args.sample_number )
        track_names = default_track_names

    return update_trace_fname, track_names

# TODO: rewrite using utils.py
def main():
    args = parse_args()

    smo = StatmapOutput( args.output_directory,
                         load_mapped_reads=True,
                         load_nc=args.use_nc )

    agg_fn = get_aggregate_fn( args.aggregate_type )
    update_trace_fname, track_names = set_track_names( smo, args )

    # specify data structures to work with
    if args.use_nc:
        mpd_rdb = smo.NC_mpd_rdb
        cond_prbs_db = smo.NC_cond_prbs_db
    else:
        mpd_rdb = smo.mpd_rdb
        cond_prbs_db = smo.cond_prbs_db

    # build the first trace from sample_number; use as the trace to update
    update_trace = load_c_trace_from_file( update_trace_fname )
    update_cond_prbs_from_trace_and_assay_type(
        mpd_rdb, cond_prbs_db, update_trace, smo.genome,
        smo.config.contents.assay_type )

    # build aggregate two traces at a time
    print >> sys.stderr, "Bootstrapping %i samples" % args.num_bootstrap_samples
    curr_trace = init_trace( smo.genome, 2, track_names )

    for i in xrange(1, args.num_bootstrap_samples):
        # NOTE: curr_trace is reset (with zero_trace) in this function
        bootstrap_trace( curr_trace, mpd_rdb, cond_prbs_db, smo.trace_update_fn )
        # aggregate bootstrapped trace into update_trace
        aggregate_over_trace_pairs( update_trace, curr_trace, agg_fn )

    close_c_traces( curr_trace )
    write_c_trace_to_file( update_trace, args.output_filename )
    close_c_traces( update_trace )

if __name__ == "__main__": main()
