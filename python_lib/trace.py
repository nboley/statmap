import sys

from ctypes import *

# load the trace library
trace_o = cdll.LoadLibrary("../src/libtrace.so")
print trace_o

print trace_o.load_trace_from_file

def load_trace( fname  ):
    print pointer( c_char_p(fname) )
    data = trace_o.load_trace_from_file( pointer( c_char_p(fname) ) )


load_trace( sys.argv[1] )
