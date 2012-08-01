import sys
import os

from ctypes import *
statmap_o = cdll.LoadLibrary( os.path.normpath( sys.path[0] +
                              "/../src/libstatmap.so" ) )

from enums import *

class c_fragment_length_dst_t(Structure):
    """
struct fragment_length_dist_t {
    int min_fl;
    int max_fl;
    float* density;
    float* chipseq_bs_density;
    float* rev_chipseq_bs_density;

    /* to user the vector operations, we need to 
       make sure that the chipseq bs densities
       are correctly aligned at 16 bytes. However, 
       we still need to free the memory so, below,
       we store the start of the allcd addresses */
    float* chipseq_bs_density_start;
    float* rev_chipseq_bs_density_start;
};
    """
    _fields_ = [
        ("min_fl", c_int),
        ("max_fl", c_int),
        ("density", POINTER(c_float)),
        ("chipseq_bs_density", POINTER(c_float)),
        ("rev_chipseq_bs_density", POINTER(c_float)),

        ("chipseq_bs_density_start", POINTER(c_float)),
        ("rev_chipseq_bs_density_start", POINTER(c_float)),
    ]

MPD_RD_ID_T = c_int

class c_mapped_reads_db_t(Structure):
    """
struct mapped_reads_db {
    FILE* fp;

    char mode; // 'r' or 'w'
    MPD_RD_ID_T num_mapped_reads;

    pthread_mutex_t* mutex;

    /* mmap data */
    /* pointer to the mmapped data and its size in bytes */
    char* mmapped_data;
    size_t mmapped_data_size; 
    
    /* mmap index */
    struct mapped_reads_db_index_t* index;
    // MPD_RD_ID_T num_indexed_reads;

    /* store the number of times that each read has been
       iterated over, and found to be below the update threshold */
    char* num_succ_iterations;
    
    struct fragment_length_dist_t* fl_dist;

    MPD_RD_ID_T current_read;
};
    """
    _fields_ = [
        ("fp", c_void_p), # unused

        ("mode", c_char),
        ("num_mapped_reads", MPD_RD_ID_T),

        ("access_lock", c_void_p), # unused (will this work?)

        ("mmapped_data", c_char_p),
        ("mmapped_data_size", c_size_t),

        ("index", c_void_p),
        
        ("num_succ_iterations", c_char_p),

        ("fl_dist", POINTER(c_fragment_length_dst_t)),

        ("current_read", MPD_RD_ID_T),
    ]

class c_cond_prbs_db_t(Structure):
    """
struct cond_prbs_db_t {
    /* this should be the same length as the number of reads in the mapped read db */
    /* each pointer should contain num_mapped_location entries */
    ML_PRB_TYPE** cond_read_prbs;
    MPD_RD_ID_T max_rd_id;
};
    """
    # note - as of Thu May 10 17:04:57 PDT 2012,
    # ML_PRB_TYPE is float
    # MPD_RD_ID_T is unsigned int
    # If either of these changes, you will need to update definitions here
    _fields_ = [
        ("cond_read_prbs", POINTER(POINTER(c_float))),
        ("max_rd_id", c_uint),
    ]

def open_mapped_reads_db_for_reading( fname ):
    c_mpd_rd_db_p = c_void_p()
    statmap_o.open_mapped_reads_db_for_reading(
            byref(c_mpd_rd_db_p), c_char_p(fname) )
    c_mpd_rd_db_p = cast( c_mpd_rd_db_p, POINTER(c_mapped_reads_db_t) )
    return c_mpd_rd_db_p

def open_mapped_reads_db_for_writing( fname ):
    c_mpd_rd_db_p = c_void_p()
    statmap_o.open_mapped_reads_db_for_writing(
            byref(c_mpd_rd_db_p), c_char_p(fname) )
    c_mpd_rd_db_p = cast( c_mpd_rd_db_p, POINTER(c_mapped_reads_db_t) )
    return c_mpd_rd_db_p

def open_mapped_reads_db_for_writing( fname ):
    c_mpd_rd

def build_fl_dist_from_filename( mpd_rd_db_p, filename ):
    statmap_o.build_fl_dist_from_filename( mpd_rd_db_p, filename )

def build_chipseq_bs_density( fl_dist_p ):
    statmap_o.build_chipseq_bs_density( fl_dist_p )

def init_cond_prbs_db_from_mpd_rdb( mpd_rd_db ):
    '''
    Return a pointer to c_cond_prbs_db_t initialized from mpd_rd_db
    '''
    c_cond_prbs_db_p = c_void_p()
    statmap_o.init_cond_prbs_db_from_mpd_rdb(
            byref(c_cond_prbs_db_p),
            mpd_rd_db )
    c_cond_prbs_db_p = cast( c_cond_prbs_db_p, POINTER(c_cond_prbs_db_t) )
    return c_cond_prbs_db_p

def reset_all_read_cond_probs( mpd_rdb, cond_prbs_db ):
    statmap_o.reset_all_read_cond_probs( mpd_rdb, cond_prbs_db )

def test():
    '''
    Test the functions in this file
    '''
    mpd_rd_db_p = open_mapped_reads_db( sys.argv[1] )
    mmap_mapped_reads_db( mpd_rd_db_p )
    index_mapped_reads_db( mpd_rd_db_p )
    cp_db_p = init_cond_prbs_db_from_mpd_rdb( mpd_rd_db_p )
    reset_all_read_cond_probs( mpd_rd_db_p, cp_db_p )

if __name__ == "__main__": test()
