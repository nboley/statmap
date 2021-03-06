import sys
import os

from ctypes import *
statmap_o = cdll.LoadLibrary( os.path.normpath( sys.path[0] +
                              "/../src/libstatmap.so" ) )

from enums import *

class c_diploid_map_data_t(Structure):
    """
struct diploid_map_data_t {
    char* chr_name;
    unsigned int chr_lens[2];
    size_t num_mappings;
    struct diploid_mapping_t* mappings;
    size_t index_len;
    struct loc_and_index_t* index;
};
    """
    _fields_ = [
        ("chr_name", c_char_p),
        ("chr_lens", c_uint*2),
        ("num_mappings", c_size_t),
        ("mappings", c_void_p), # unused
        ("index_len", c_size_t),
        ("index", c_void_p), # unused
    ]

class c_diploid_maps_t(Structure):
    """
struct diploid_maps_t {
    int count;
    struct diploid_map_data_t* maps;
};

    """
    _fields_ = [
        ("count", c_int),
        ("maps", POINTER(c_diploid_map_data_t)),
    ]

class c_index_t(Structure):
    """
struct index_t {
    void* index;
    enum INDEX_TYPE index_type;
    /* the length of a sequence in the index */
    int seq_length;

    /* the pseudo locations info */
    struct pseudo_locations_t* ps_locs;
    
    /* diploid map data */
    int num_diploid_chrs;
    /* array of dipoid_map_data_t */
    struct diploid_map_data_t* map_data;
};
    """
    _fields_ = [
        ("index", c_void_p),
        ("index_type", ENUM_TYPE),
        ("seq_length", c_int),

        ("ps_locs", c_void_p),

        ("diploid_maps", POINTER(c_diploid_maps_t)),
    ]

class c_genome_data_t(Structure):
    """
struct genome_data {
    /* The index for this genome */
    struct index_t* index;

    /* The number of chromosomes */
    int num_chrs;
    /* The chromosome names - indexed by their keys */
    char** chr_names;
    /* chromosomes - the actual sequences */
    char** chrs;
    /* the length of the chrsomosomes in bps' */
    unsigned int* chr_lens;
    /* the source the chromosomes */
    enum CHR_SOURCE* chr_sources;

    enum bool is_mmapped;
};
    """
    _fields_ = [
        ("index", POINTER(c_index_t)),

        ("num_chrs", c_int),
        ("chr_names", POINTER(c_char_p)),
        ("chrs", POINTER(c_char_p)),
        ("chr_lens", POINTER(c_int)),
        ("chr_sources", POINTER(ENUM_TYPE)), # enum CHR_SOURCE

        ("is_mmapped", ENUM_TYPE), # enum bool
    ]

### Functions ###

def init_genome():
    c_genome_p = c_void_p()
    statmap_o.init_genome( byref(c_genome_p) )
    c_genome_p = cast( c_genome_p, POINTER(c_genome_data_t ) )
    return c_genome_p

def load_genome_from_disk( fname ):
    '''
    Given the name of a genome file, returns a pointer to c_genome_data_t
    '''
    c_genome_data_p = c_void_p()
    statmap_o.load_genome_from_disk( byref(c_genome_data_p), c_char_p(fname) )
    c_genome_data_p = cast( c_genome_data_p, POINTER(c_genome_data_t) )
    return c_genome_data_p

def add_chrs_from_fasta_file( genome, filename, chr_source ):
    statmap_o.add_chrs_from_fasta_file( genome, filename, chr_source )

def init_index( index, seq_len, num_diploid_chrs ):
    statmap_o.init_index( index, seq_len, num_diploid_chrs )

def index_genome( genome ):
    statmap_o.index_genome( genome )

def write_genome_to_disk( genome, output_fname ):
    statmap_o.write_genome_to_disk( genome, output_fname )

def build_ondisk_index( index, index_fname ):
    statmap_o.build_ondisk_index( index, index_fname )

def test():
    '''
    Test the functions in this file
    '''
    # path to a genome file
    genome = load_genome_from_disk( sys.argv[1] )

    # print the # of chrs
    print genome.contents.num_chrs

    # print the names of all the chrs
    for chr in range(genome.contents.num_chrs):
        print genome.contents.chr_names[chr]

if __name__ == "__main__": test()
