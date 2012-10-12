import sys
import os

from ctypes import *
statmap_o = cdll.LoadLibrary( os.path.normpath( sys.path[0] +
                              "/../src/libstatmap.so" ) )

from enums import *

class c_args_t(Structure):
    """
struct args_t {
    char* genome_fname;
    char* genome_index_fname;
    
    char* unpaired_reads_fnames;
    char* pair1_reads_fnames;
    char* pair2_reads_fnames;
    struct rawread_db_t* rdb;

    char* unpaired_NC_reads_fnames;
    char* pair1_NC_reads_fnames;
    char* pair2_NC_reads_fnames;
    struct rawread_db_t* NC_rdb;

    char* frag_len_fname;
    FILE* frag_len_fp;

    char* output_directory;

    char* sam_output_fname;

    float mapping_metaparameter;

    int min_num_hq_bps;

    int num_starting_locations;

    int num_threads;

    enum error_model_type_t error_model_type;
        
    enum input_file_type_t input_file_type;
    enum assay_type_t assay_type; 

    int max_reference_insert_len;
    int softclip_len;
};
    """
    _fields_ = [
        ("genome_fname", c_char_p),
        ("genome_index_fname", c_char_p),

        ("unpaired_reads_fnames", c_char_p),
        ("pair1_reads_fnames", c_char_p),
        ("pair2_reads_fnames", c_char_p),
        ("rdb", c_void_p),

        ("unpaired_NC_reads_fnames", c_char_p),
        ("pair1_NC_reads_fnames", c_char_p),
        ("pair2_NC_reads_fnames", c_char_p),
        ("NC_rdb", c_void_p),

        ("frag_len_fname", c_char_p),
        ("frag_len_fp", c_void_p), # really, why would we save a file* if we have the fname?

        ("output_directory", c_char_p),

        ("sam_output_fname", c_char_p),
        
        ("mapping_metaparameter", c_float),

        ("min_num_hq_bps", c_int),

        ("num_starting_locations", c_int),

        ("num_threads", c_int),

        ("error_model_type", ENUM_TYPE),   # enum error_model_type_t
        
        ("input_file_type", ENUM_TYPE),    # enum input_file_type_t
        ("assay_type", ENUM_TYPE),         # enum assay_type_t

        ("max_reference_insert_len", c_int),
        ("softclip_len", c_int),
    ]

def load_config_from_file( fname ):
    '''
    Given the fname of a config.dat file, returns a pointer to c_args_t
    '''
    c_args_p = c_void_p()
    statmap_o.read_config_file_fname_from_disk( c_char_p(fname), byref(c_args_p) )
    c_args_p = cast( c_args_p, POINTER(c_args_t) )
    return c_args_p

def test():
    '''
    Test the functions in this file
    '''
    args = load_config_from_file( sys.argv[1] )
    print args.contents.genome_fname

    # test working with global variables
    print statmap_o.get_num_threads()
    print statmap_o.get_min_num_hq_bps()

    # set global variables
    statmap_o.set_num_threads( 4 )
    statmap_o.set_min_num_hq_bps( 10 )

    # print them again
    print statmap_o.get_num_threads()
    print statmap_o.get_min_num_hq_bps()

if __name__ == "__main__": test()
