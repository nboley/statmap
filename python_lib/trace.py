import sys

import numpy

from collections import OrderedDict

from ctypes import *
# load the trace library
trace_o = cdll.LoadLibrary("../src/libtrace.so")

class c_trace_t(Structure):
    """
#define TRACE_TYPE float
    
struct trace_t {
    int num_tracks;
    char** track_names;
    
    int num_chrs;
    char** chr_names;
    unsigned int* chr_lengths;
    
    /* num_traces X num_chrs matrix */
    TRACE_TYPE*** traces;

    #ifdef USE_MUTEX
    pthread_mutex_t*** locks;
    #else
    pthread_spinlock_t*** locks;
    #endif
};
"""
    _fields_ = [("num_tracks", c_int),
                ("track_names", POINTER(c_char_p)),

                ("num_chrs", c_int),
                ("chr_names", POINTER(c_char_p)),
                ("chr_lengths", POINTER(c_uint)),
                
                ("traces", POINTER(POINTER(POINTER(c_float)))),
                
                # we should never touch this from python, 
                # but I need it to make sizes correct 
                ("locks", c_void_p )
               ]

class trace_t(OrderedDict):
    def __init__( self, c_trace_p ):
        # call the parent's initializer
        OrderedDict.__init__(self)
        
        # de-reference the trace
        c_trace = c_trace_p[0]
        
        for track_index in range(c_trace.num_tracks):
            track_name = c_trace.track_names[track_index]
            self[track_name] = OrderedDict()
            for chr_index in range(c_trace.num_chrs):
                chr_name, chr_len = c_trace.chr_names[chr_index], c_trace.chr_lengths[chr_index]
                self[track_name][chr_name] = numpy.zeros( chr_len  )
                # Make this faster - I must be able to load this through a pointer copy
                for index in range( c_trace.chr_lengths[chr_index] ):
                    self[track_name][chr_name][index] \
                        = c_trace.traces[track_index][chr_index][index]
        return

def load_c_trace_from_file( fname  ):
    c_trace_p = c_void_p()    
    data = trace_o.load_trace_from_file( byref(c_trace_p), fname )
    c_trace_p = cast( trace_p, POINTER(c_trace_t) )
    return c_trace_p

def load_trace_from_file( fname  ):
    c_trace_p = c_void_p()    
    data = trace_o.load_trace_from_file( byref(c_trace_p), c_char_p(fname) )
    print( data )
    c_trace_p = cast( c_trace_p, POINTER(c_trace_t) )
    
    return trace_t( c_trace_p )


trace = load_trace_from_file( sys.argv[1] )
print( trace )
