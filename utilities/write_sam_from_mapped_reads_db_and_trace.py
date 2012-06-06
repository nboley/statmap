#!/usr/bin/python

import sys
import os

# add python_lib to sys.path
sys.path.insert(0, os.path.normpath( sys.path[0] + "/../python_lib") )
from sam import *
from trace import *
from utils import *

def usage():
    print "Usage: ./mapped_reads_to_sam output_directory input.bin.trace [--use-nc]"
    sys.exit(1)

def parse_args():
    import argparse; parser = argparse.ArgumentParser()

    # Required arguments
    parser.add_argument("output_directory",
        help="Path to Statmap output directory to process")
    parser.add_argument("input_trace",
        help="Binary trace file to use to set conditional probabilities")

    # Optional arguments
    parser.add_argument("--use-nc", "-nc",
        help="(For Chipseq only) Use the negative control mapped reads.",
        action="store_true")

    return parser.parse_args()

def main():
    args = parse_args()
    it_path = os.path.abspath( args.input_trace ) # use absolute path

    smo = StatmapOutput( sys.argv[1],
                         load_mapped_reads=True,
                         load_raw_reads=True,
                         load_nc=args.use_nc )

    # specify data structures to work with (using nc vs ip/normal)
    if args.use_nc:
        raw_rdb = smo.NC_rawread_db
        mpd_rdb = smo.NC_mpd_rdb
        cond_prbs_db = smo.NC_cond_prbs_db
    else:
        raw_rdb = smo.rawread_db
        mpd_rdb = smo.mpd_rdb
        cond_prbs_db = smo.cond_prbs_db

    # default not to reset conditional probabilities
    reset_cond_prbs = False

    trace = load_c_trace_from_file( it_path )

    # update the read conditional probabilities
    update_cond_prbs_from_trace_and_assay_type(
        mpd_rdb, cond_prbs_db, trace, smo.genome,
        smo.config.contents.assay_type )

    close_c_traces( trace )

    # write mapped reads to sam
    write_mapped_reads_to_sam(
        raw_rdb, mpd_rdb, cond_prbs_db, smo.genome,
        reset_cond_prbs, False,
        sys.stdout )

if __name__ == "__main__": main()
