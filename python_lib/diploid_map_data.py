import sys
import os

from ctypes import *
statmap_o = cdll.LoadLibrary( os.path.normpath( sys.path[0] +
                              "/../src/libstatmap.so" ) )

from enums import *

def get_chr_index( genome, chr_name ):
    return statmap_o.find_chr_index( genome, chr_name )

def parse_map_file( map_fname, map_data_t, genome,
                    paternal_chr_index, maternal_chr_index ):
    statmap_o.parse_map_file( map_fname, map_data_t, genome,
            paternal_chr_index, maternal_chr_index )

def index_diploid_map_data( map_data_t ):
    statmap_o.index_diploid_map_data( map_data_t )
