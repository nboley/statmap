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

class c_mapped_reads_db_index_t(Structure):
    """
struct mapped_reads_db_index_t {
    MPD_RD_ID_T read_id;
    char* ptr;
};
    """
    _fields_ = [
        ("read_id", MPD_RD_ID_T),
        ("ptr", c_char_p),
    ]

class c_mapped_reads_db_t(Structure):
    """
struct mapped_reads_db {

    /* A read may map, be considered unmappable, or be considered mappable but
     * fail to map (nonmapping). These files store data about each type of
     * read. */
    FILE* mapped_fp;
    MPD_RD_ID_T num_mapped_reads;
    pthread_mutex_t* mapped_mutex;

    FILE* unmappable_fp;
    MPD_RD_ID_T num_unmappable_reads;
    pthread_mutex_t* unmappable_mutex;

    FILE* nonmapping_fp;
    MPD_RD_ID_T num_nonmapping_reads;
    pthread_mutex_t* nonmapping_mutex;

    char mode; // 'r' or 'w'

    /* mmap data */
    /* pointer to the mmapped data and its size in bytes */
    char* mmapped_data;
    size_t mmapped_data_size; 
    
    /* mmap index */
    struct mapped_reads_db_index_t* index;
    
    /* store the number of times that each read has been
       iterated over, and found to be below the update threshold */
    char* num_succ_iterations;
    
    struct fragment_length_dist_t* fl_dist;

    MPD_RD_ID_T current_read;
};
    """
    _fields_ = [
        ("mapped_fp", c_void_p),
        ("num_mapped_reads", MPD_RD_ID_T ),
        ("mapped_mutex", c_void_p),

        ("unmappable_fp", c_void_p),
        ("num_unmappable_reads", MPD_RD_ID_T ),
        ("unmappable_mutex", c_void_p),

        ("nonmapping_fp", c_void_p),
        ("num_nonmapping_reads", MPD_RD_ID_T ),
        ("nonmapping_mutex", c_void_p),

        ("mode", c_char),

        ("mmapped_data", c_char_p),
        ("mmapped_data_size", c_size_t),

        ("index", POINTER(c_mapped_reads_db_index_t)),
        
        ("num_succ_iterations", c_char_p),

        ("fl_dist", POINTER(c_fragment_length_dst_t)),

        ("current_read", MPD_RD_ID_T),
    ]

ML_PRB_TYPE = c_float

class c_cond_prbs_db_t(Structure):
    """
struct cond_prbs_db_t {
    /* this should be the same length as the number of reads in the mapped read db */
    /* each pointer should contain num_mapped_location entries */
    ML_PRB_TYPE** cond_read_prbs;
    MPD_RD_ID_T max_rd_id;
};
    """
    _fields_ = [
        ("cond_read_prbs", POINTER(POINTER(ML_PRB_TYPE))),
        ("max_rd_id", MPD_RD_ID_T),
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
    mpd_rd_db_p = open_mapped_reads_db_for_reading( sys.argv[1] )
    cp_db_p = init_cond_prbs_db_from_mpd_rdb( mpd_rd_db_p )
    reset_all_read_cond_probs( mpd_rd_db_p, cp_db_p )

if __name__ == "__main__": test()
