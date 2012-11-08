import sys
import os

from ctypes import *
statmap_o = cdll.LoadLibrary( os.path.normpath( sys.path[0] +
                              "/../src/libstatmap.so" ) )

from enums import *

def init_logging():
    statmap_o.init_logging()

# FIXME - figure out how to call this on close for StatmapOutput
# For now, since we're just using stderr, it will be handled automatically by the OS
def finish_logging():
    statmap_o.finish_logging()
