import sys

from ctypes import *

# load the trace library
trace_o = cdll.LoadLibrary("../src/libtrace.so")


def load_trace( fname  ):
    trace = c_void_p()
    data = trace_o.load_trace_from_file( byref(trace), fname )


load_trace( sys.argv[1] )
