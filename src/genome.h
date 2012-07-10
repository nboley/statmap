/* Copyright (c) 2009-2010 Nathan Boley */

#ifndef GENOME
#define GENOME

#include "config.h"

/* 
 * Eventually we may want indexes to pluggable, providing for multiple types ( or
 * for instance, some on disk, etc. ). To this end, we define some generic index
 * infrastructure.
 *
 */

enum INDEX_TYPE {
    TREE = 1
};

struct index_t {
    void* index;
    enum INDEX_TYPE index_type;
    /* the length of a sequence in the index */
    int seq_length;

    /* the pseudo locations info */
    struct pseudo_locations_t* ps_locs;

    /* diploid map data */
    struct diploid_maps_t* diploid_maps;    
};

struct static_node;

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

void
add_chr_to_genome( char* chr_name, char* chr_str, unsigned int chr_len, enum CHR_SOURCE chr_source, struct genome_data* gen ); 

int
find_chr_index( struct genome_data* genome, const char* const chr_name );

int
find_diploid_chr_index(
    struct genome_data* genome,
    const char* const prefix,
    enum CHR_SOURCE source
);

char* get_chr_prefix( char* chr_name );

int get_map_data_index_from_chr_index(
    struct genome_data* genome,
    int chr_index
);

char* 
find_seq_ptr( struct genome_data* genome, 
              int chr_index, unsigned int loc, 
              int read_len );

extern void
index_genome( struct genome_data* genome );

extern void 
add_chrs_from_fasta_file( struct genome_data* gen, char* fname, enum CHR_SOURCE chr_source );

void
init_genome( struct genome_data** gen );

void
free_genome( struct genome_data* gen );

struct genome_header {
    size_t size;
    size_t genome_offset;
    size_t pseudo_locs_offset;
    size_t index_offset;
};

void
load_genome_from_disk( struct genome_data** gen, char* fname );

void
write_genome_to_disk( struct genome_data* gen, char* fname );


#endif //#ifndef GENOME
