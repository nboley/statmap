import sys

from ctypes import *
statmap_o = cdll.LoadLibrary("../src/libstatmap.so")

from enums import *

def call_peaks( genome ):
    statmap_o.call_peaks( genome )

def call_peaks_at_local_maxima( genome, samples_dname ):
    statmap_o.call_peaks_at_local_maxima( genome, samples_dname )
