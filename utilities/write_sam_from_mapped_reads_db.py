# assume this is in utilities
import sys
import os

# add python_lib to sys.path
sys.path.insert(0, os.path.normpath( sys.path[0] + "/../python_lib") )
from sam import *
from trace import *
from utils import *

def parse_args():
    import argparse; parser = argparse.ArgumentParser()

    # Required arguments
    parser.add_argument("output_directory",
        help="Path to Statmap output directory to process")

    # Optional arguments
    parser.add_argument("--sample-num", "-s",
        help="Sample number to get conditional probabilities from",
        type=int,
        default=0)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--use-ip", "-ip",
        help="(For Chipseq only) Use the immunoprecipitate mapped reads (default). ",
        action="store_true")
    group.add_argument("--use-nc", "-nc",
        help="(For Chipseq only) Use the negative control mapped reads. ",
        action="store_true")

    return parser.parse_args()

def main():
    args = parse_args()

    smo = StatmapOutput( sys.argv[1],
                         load_mapped_reads=True,
                         load_raw_reads=True,
                         load_nc=args.use_nc )

    # if Chip-Seq, require the user to specify which reads to use
    if smo.config.contents.assay_type == CHIP_SEQ and not \
            (args.use_ip or args.use_nc):
        print >> sys.stderr, "ERROR    : No Chip-seq reads specified - expecting --use-nc ( for use negative control ) or --use-ip ( for immunoprecipitate )"
        print >> sys.stderr, "HINT     : This experiment used a negative control - you need to specify which reads to output."
        sys.exit(1)

    # default to reset conditional probabilities
    reset_cond_prbs = True

    # specify data structures to work with
    if args.use_nc:
        raw_rdb = smo.NC_rawread_db
        mpd_rdb = smo.NC_mpd_rdb
        cond_prbs_db = smo.NC_cond_prbs_db
    else:
        raw_rdb = smo.rawread_db
        mpd_rdb = smo.mpd_rdb
        cond_prbs_db = smo.cond_prbs_db

    # if we specified a marginal density
    if args.sample_num > 0:
        reset_cond_prbs = False

        # load the trace that stores the marginal read density we're
        # interested in
        if not (args.use_ip or args.use_nc):
            trace_fname = "./samples/sample%i.bin.trace" % args.sample_num
        elif args.use_ip:
            trace_fname = "./samples/sample%i.ip.bin.trace" % args.sample_num
        elif args.use_nc:
            trace_fname = "./samples/sample%i.nc.bin.trace" % args.sample_num

        trace = load_c_trace_from_file( trace_fname )
        # update the read conditional probabilities based on the assay
        # and the selected trace
        update_cond_prbs_from_trace_and_assay_type(
            mpd_rdb, cond_prbs_db, trace, smo.genome,
            smo.config.contents.assay_type
        )

    # write mapped reads to sam
    write_mapped_reads_to_sam(
        raw_rdb, mpd_rdb, cond_prbs_db, smo.genome,
        reset_cond_prbs, False,
        sys.stdout
    )
         
if __name__ == "__main__": main()
