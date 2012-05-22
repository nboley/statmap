import sys

from ctypes import *
statmap_o = cdll.LoadLibrary("../src/libstatmap.so")

from enums import *

class c_genome_data_t(Structure):
    """
struct genome_data {
    /* The index for this genome */
    struct index_t* index;

    /* the pseudo locations info */
    struct pseudo_locations_t* ps_locs;
    
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
        ("index", c_void_p), # unused

        ("ps_locs", c_void_p), # unused

        ("num_chrs", c_int),
        ("chr_names", POINTER(c_char_p)),
        ("chrs", POINTER(c_char_p)),
        ("chr_lens", POINTER(c_int)),
        ("chr_sources", POINTER(c_uint)), # enum CHR_SOURCE

        ("is_mmapped", c_uint), # enum bool
    ]

### Functions ###

def load_genome_from_disk( fname ):
    '''
    Given the name of a genome file, returns a pointer to c_genome_data_t
    '''
    c_genome_data_p = c_void_p()
    statmap_o.load_genome_from_disk( byref(c_genome_data_p), c_char_p(fname) )
    c_genome_data_p = cast( c_genome_data_p, POINTER(c_genome_data_t) )
    return c_genome_data_p

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
