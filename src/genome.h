/* Copyright (c) 2009-2010 Nathan Boley */

#ifndef GENOME
#define GENOME

/* 
 * Eventually we may want indexes to pluggable, providing for multiple types ( or
 * for instance, some on disk, etc. ). To this end, we define some generic index
 * infrastructure.
 *
 */

enum INDEX_TYPE {
    TREE = 1
};

typedef struct {
    void* index;
    enum INDEX_TYPE index_type;
    /* the length of a sequence in the index */
    int seq_length;
} index_t;

struct static_node;
struct snp_db_t;

struct genome_data {
    /* The index for this genome */
    index_t* index;
    
    /* snps, NULL if there are none */
    struct snp_db_t* snp_db;
    
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
};

void
add_chr_to_genome( char* chr_name, char* chr_str, unsigned int chr_len,  struct genome_data* gen ); 

int
find_chr_index( struct genome_data* genome, const char* const chr_name );

extern void 
add_chrs_from_fasta_file( struct genome_data* gen, FILE* f  );

void
init_genome( struct genome_data** gen );

void
free_genome( struct genome_data* gen );


#endif //#ifndef GENOME
