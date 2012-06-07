import sys
import os

from ctypes import *
statmap_o = cdll.LoadLibrary( os.path.normpath( sys.path[0] +
                              "/../src/libstatmap.so" ) )

from enums import *

class c_rawread_db_t(Structure):
    """
struct rawread_db_t {
    /* 
       This 'db' is functioning as both a cursor and a db. Therefore, 
       readkey stores the key of the current read. Calling get_next_read
       will pull out the next read, and increment this. rewind will
       reset this to 0. 
    */
    readkey_t readkey;
    /* We lock this mutex whenever we grab a read */
    pthread_spinlock_t* lock; 

    FILE* single_end_reads;
    FILE* paired_end_1_reads;
    FILE* paired_end_2_reads;

    FILE* unmappable_single_end_reads;
    FILE* unmappable_paired_end_1_reads;
    FILE* unmappable_paired_end_2_reads;

    FILE* non_mapping_single_end_reads;
    FILE* non_mapping_paired_end_1_reads;
    FILE* non_mapping_paired_end_2_reads;

    enum inputfile_type file_type;
    enum assay_type_t assay;
};
    """
    _fields_ = [
        ("readkey", c_uint), # readkey_t is unsigned int
        
        ("lock", c_void_p), # unused - pthread_spinlock_t*

        ("single_end_reads", c_void_p),
        ("paired_end_1_reads", c_void_p),
        ("paired_end_2_reads", c_void_p),
        
        ("unmappable_single_end_reads", c_void_p),
        ("unmappable_paired_end_1_reads", c_void_p),
        ("unmappable_paired_end_2_reads", c_void_p),

        ("non_mapping_single_end_reads", c_void_p),
        ("non_mapping_paired_end_1_reads", c_void_p),
        ("non_mapping_paired_end_2_reads", c_void_p),

        ("file_type", c_uint),  # enum inputfile_type
        ("assay", c_uint),      # enum assay_type_t
    ]

def populate_rawread_db( unpaired_fname,
                         pair1_fname,
                         pair2_fname ):
    """
    Return a pointer to a populated rawread db. Loads single ended or paired
    end reads, depending on what files are available.
    """
    # init rawread db
    rawread_db_p = c_void_p()
    statmap_o.init_rawread_db( byref(rawread_db_p) )

    # try loading unpaired reads
    try:
        tmp_fp = open( unpaired_fname )
        tmp_fp.close()
        statmap_o.add_single_end_reads_to_rawread_db(
            rawread_db_p, unpaired_fname, FASTQ, UNKNOWN
        )
    except:
        pass

    # try loading paired end reads
    try:
        tmp_fp = open( pair1_fname )
        tmp_fp.close()
        statmap_o.add_paired_end_reads_to_rawread_db(
            rawread_db_p, pair1_fname, pair2_fname, FASTQ, UNKNOWN
        )
    except:
        pass

    rawread_db_p = cast( rawread_db_p, POINTER(c_rawread_db_t) )
    return rawread_db_p

def test():
    """
    test functions in this file
    """
    pass

if __name__ == "__main__": test()
