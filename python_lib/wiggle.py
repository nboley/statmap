import sys
import os

from ctypes import *
statmap_o = cdll.LoadLibrary( os.path.normpath( sys.path[0] +
                              "/../src/libstatmap.so" ) )

from enums import *
from utils import *

def aggregate_over_wiggles( wigs, num_wigs, ofp, threshold, agg_fn ):
    # convert list of Python File objects wigs -> array of FILE*
    num_wigs_c_array = FILE_ptr * num_wigs
    wigs_c_array = num_wigs_c_array( *[ PyFile_AsFile(fobj) for fobj in wigs ] )

    statmap_o.aggregate_over_wiggles(
        wigs_c_array, num_wigs,
        PyFile_AsFile( ofp ),
        c_float( threshold ),
        agg_fn )

def write_wiggle_from_trace_to_stream( trace, os, filter_threshold ):
    # converts Python File object os -> FILE*
    statmap_o.write_wiggle_from_trace_to_stream( trace,
        PyFile_AsFile(os),
        c_double( filter_threshold ) )
