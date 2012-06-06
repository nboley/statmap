import sys
import os

# add python_lib to sys.path (at the front, since there's a trace module in 
# the standard library )
sys.path.insert(0, "../python_lib/" )

from config_parsing import *
from genome import *
from mapped_read import *
from trace import *
from utils import *

### CONSTANTS - from config.h ###
MAX_NUM_EM_ITERATIONS           = 500
MAX_PRB_CHANGE_FOR_CONVERGENCE  = 1e-2
BOOTSTRAP_SAMPLES_ALL_PATH      = "./bootstrap_samples/all_traces/"

STARTING_SAMPLES_META_INFO_FNAME = "./starting_samples/meta_info.csv"
RELAXED_SAMPLES_META_INFO_FNAME  = "./samples/meta_info.csv"

def usage():
    print "Usage: ./sample_mapping.py output_directory sample_number [num_threads]"
    sys.exit(1)

def determine_next_sample_index():
    # for now, determine from command line argument
    return int( sys.argv[2] )

def main():
    if len(sys.argv) not in (3,4): usage()

    if len(sys.argv) == 4:  # if the number of threads to use is given
        num_threads = int( sys.argv[3] )
    else:                   # otherwise let StatmapOutput determine it
        num_threads = 0

    smo = StatmapOutput( sys.argv[1], num_threads,
                         load_mapped_reads=True,
                         load_raw_reads=True,
                         load_nc=True )

    # open the meta data file to append to
    meta_info_fp = open( RELAXED_SAMPLES_META_INFO_FNAME, "a" )
    ss_meta_info_fp = open( STARTING_SAMPLES_META_INFO_FNAME, "a" )

    # determine sample index
    sample_index = determine_next_sample_index()

    use_random_start = True
    max_prb_change_for_convergence = MAX_PRB_CHANGE_FOR_CONVERGENCE

    if smo.config.contents.assay_type == CHIP_SEQ:
        if smo.NC_mpd_rdb:
            statmap_o.take_chipseq_sample_wnc(
                smo.mpd_rdb, smo.NC_mpd_rdb,
                smo.genome,
                PyFile_AsFile(meta_info_fp), # -> FILE*
                sample_index,
                c_float(max_prb_change_for_convergence),
                use_random_start
            )
        else:
            statmap_o.take_chipseq_sample(
                smo.mpd_rdb,
                smo.genome,
                PyFile_AsFile(ss_meta_info_fp), # -> FILE*
                PyFile_AsFile(meta_info_fp), # -> FILE*
                MAX_NUM_EM_ITERATIONS,
                c_float(max_prb_change_for_convergence)
            )
    elif smo.config.contents.assay_type == CAGE:
        statmap_o.take_cage_sample(
            smo.mpd_rdb,
            smo.genome,
            sample_index,
            PyFile_AsFile(ss_meta_info_fp), # -> FILE*
            PyFile_AsFile(meta_info_fp), # -> FILE*
            MAX_NUM_EM_ITERATIONS,
            c_float(max_prb_change_for_convergence)
        )
    else:
        print >> sys.stderr, "FATAL    : Unsupported assay type %i" \
                % smo.config.contents.assay_type

if __name__ == "__main__": main()
