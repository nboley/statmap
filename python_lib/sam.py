import sys

import ctypes
statmap_o = ctypes.cdll.LoadLibrary("../src/libstatmap.so")

from enums import *
from utils import *

def write_mapped_reads_to_sam( rdb, mappings_db, cond_prbs_db, genome,
                               reset_cond_read_prbs, expand_pseudo_locations,
                               sam_ofp ):
    statmap_o.write_mapped_reads_to_sam(
        rdb, mappings_db, cond_prbs_db, genome,
        reset_cond_read_prbs, expand_pseudo_locations,
        # convert Python file object to C FILE*
        PyFile_AsFile( sam_ofp )
    )

def test():
    smo = StatmapOutput( sys.argv[1],
                         load_mapped_reads=True,
                         load_rawreads=True )
    write_mapped_reads_to_sam(
        smo.rawread_db, smo.mpd_rdb, smo.cond_prbs_db, smo.genome,
        True, False,
        sys.stdout
    )

if __name__ == "__main__": test()
