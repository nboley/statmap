import sys
import os

from ctypes import *
statmap_o = cdll.LoadLibrary( os.path.normpath( sys.path[0] +
                              "/../src/libstatmap.so" ) )

import numpy

try:
    from collections import OrderedDict
except ImportError:
    from collections import MutableMapping

    class OrderedDict(dict, MutableMapping):
        # Methods with direct access to underlying attributes

        def __init__(self, *args, **kwds):
            if len(args) > 1:
                raise TypeError('expected at 1 argument, got %d', len(args))
            if not hasattr(self, '_keys'):
                self._keys = []
            self.update(*args, **kwds)

        def clear(self):
            del self._keys[:]
            dict.clear(self)

        def __setitem__(self, key, value):
            if key not in self:
                self._keys.append(key)
            dict.__setitem__(self, key, value)

        def __delitem__(self, key):
            dict.__delitem__(self, key)
            self._keys.remove(key)

        def __iter__(self):
            return iter(self._keys)

        def __reversed__(self):
            return reversed(self._keys)

        def popitem(self):
            if not self:
                raise KeyError
            key = self._keys.pop()
            value = dict.pop(self, key)
            return key, value

        def __reduce__(self):
            items = [[k, self[k]] for k in self]
            inst_dict = vars(self).copy()
            inst_dict.pop('_keys', None)
            return (self.__class__, (items,), inst_dict)

        # Methods with indirect access via the above methods

        setdefault = MutableMapping.setdefault
        update = MutableMapping.update
        pop = MutableMapping.pop
        keys = MutableMapping.keys
        values = MutableMapping.values
        items = MutableMapping.items

        def __repr__(self):
            pairs = ', '.join(map('%r: %r'.__mod__, self.items()))
            return '%s({%s})' % (self.__class__.__name__, pairs)

        def copy(self):
            return self.__class__(self)

        @classmethod
        def fromkeys(cls, iterable, value=None):
            d = cls()
            for key in iterable:
                d[key] = value
            return d

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
        
        # keep a reference to c_trace_p?

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

def init_trace( genome_p, num_tracks, track_names ):
    '''
    Init a trace and return it
    '''

    # convert python list of strings into C array of char*
    c_track_names = (c_char_p * num_tracks)()
    for i, c_track in enumerate(c_track_names):
        c_track_names[i] = c_char_p( track_names[i] )

    c_trace_p = c_void_p()
    statmap_o.init_trace( genome_p, byref(c_trace_p), num_tracks, c_track_names )
    c_trace_p = cast( c_trace_p, POINTER(c_trace_t) )
    return c_trace_p

def load_c_trace_from_file( fname ):
    c_trace_p = c_void_p()    
    statmap_o.load_trace_from_file( byref(c_trace_p), fname )
    c_trace_p = cast( c_trace_p, POINTER(c_trace_t) )
    return c_trace_p

def load_trace_from_file( fname ):
    c_trace_p = c_void_p()    
    statmap_o.load_trace_from_file( byref(c_trace_p), c_char_p(fname) )
    c_trace_p = cast( c_trace_p, POINTER(c_trace_t) )
    return trace_t( c_trace_p )

def write_c_trace_to_file( c_trace_p, fname ):
    '''
    Writes c_trace_t c_trace to fname using C function write_trace_to_file
    '''
    statmap_o.write_trace_to_file( c_trace_p, c_char_p(fname) )

def close_c_traces( c_trace_p ):
    '''
    Closes the trace c_trace
    '''
    statmap_o.close_traces( c_trace_p )

def update_traces_from_mapped_reads( mpd_rdb, cond_prbs_db, trace, update_fn ):
    statmap_o.update_traces_from_mapped_reads(
        mpd_rdb, cond_prbs_db, trace, update_fn
    )

def update_cond_prbs_from_trace_and_assay_type( mpd_rd_db_p, cond_prbs_db_p,
        c_trace_p, genome_p, assay_type ):
    statmap_o.update_cond_prbs_from_trace_and_assay_type(
        mpd_rd_db_p, cond_prbs_db_p, c_trace_p, genome_p, assay_type
    )

def bootstrap_trace( c_trace_p,
                     mpd_rd_db_p,
                     cond_prbs_db_p,
                     update_trace_expectation_from_location ):
    statmap_o.bootstrap_traces_from_mapped_reads(
            mpd_rd_db_p,
            cond_prbs_db_p,
            c_trace_p,
            update_trace_expectation_from_location )

def aggregate_over_trace_pairs( update_trace_p, curr_trace_p, agg_fn ):
    statmap_o.aggregate_over_trace_pairs( update_trace_p, curr_trace_p, agg_fn )

if __name__ == "__main__":
    trace = load_trace_from_file( sys.argv[1] )
    print( trace )
