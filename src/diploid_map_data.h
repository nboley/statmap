#ifndef DIPLOID_MAP
#define DIPLOID_MAP

#include "config.h"

/*
 * Stores mapping from a single line in .map file
 */
struct diploid_mapping_t {
    SIGNED_LOC paternal;
    SIGNED_LOC maternal;
};

/*
 * Location of a haplotype start (basepair number)
 * and the index of the corresponding entry in diploid_map_data_t->mappings
 */
struct loc_and_index_t {
    SIGNED_LOC loc;
    int index;
};

/*
 * Stores all informaiton parsed from a .map file, as well as
 * the index of loc_and_index_t types that are used to find the index
 * in the .map file of a mapping given the loc.
 */
struct diploid_map_data_t {
    char* chr_name;
    unsigned int chr_lens[2];
    size_t num_mappings;
    struct diploid_mapping_t* mappings;
    size_t index_len;
    struct loc_and_index_t* index;
};

/*
 * Represents an identical segment on both chromosomes by
 * the start of the region on the paternal chromosome,
 *      "           "           " maternal chromosme,
 * the length of the segment.
 */
struct chr_subregion_t {
    SIGNED_LOC paternal_start_pos;
    SIGNED_LOC maternal_start_pos;
    SIGNED_LOC segment_length;
};

void
init_diploid_map_data( 
    struct diploid_map_data_t** map_data, 
    char* chr_name, 
    unsigned int* chr_lens
);

void
free_diploid_map_data(
    struct diploid_map_data_t* map_data
);

void
write_diploid_map_data_to_file(
    struct diploid_map_data_t* map_data,
    int num_chrs,
    FILE* fp
);

int
read_diploid_map_data_from_file(
    struct diploid_map_data_t** map_data,
    FILE* fp
);

void
index_diploid_map_data( 
    struct diploid_map_data_t* data
);

int
find_diploid_locations( 
    struct diploid_map_data_t* data, 
    SIGNED_LOC paternal_pos
);

void
parse_map_file(
    char* fname, 
    struct diploid_map_data_t** map_data,
    struct genome_data* genome,
    int paternal_chr_index,
    int maternal_chr_index
);

void
build_unique_sequence_segments(
    struct diploid_map_data_t* data, 
    int seq_len,
    struct chr_subregion_t** segments,
    int* num_segments 
);
#endif // #ifndef DIPLOID_MAP
