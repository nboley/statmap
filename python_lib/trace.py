import sys

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
    c_trace_p = cast( c_trace_p, POINTER(c_trace_t) )
    return trace_t( c_trace_p )


if __name__ == "__main__":
    trace = load_trace_from_file( sys.argv[1] )
    print( trace )
