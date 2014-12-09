/* Copyright (c) 2009-2012 Nathan Boley */

#ifndef READ_HEADER
#define READ_HEADER

#include "config.h"
#include "rawread.h"
#include "quality.h"
#include "error_correction.h"

struct mapping_params; // fwd declaration (?)

#define POS_SINGLE_END 0
#define POS_PAIRED_END_1 0
#define POS_PAIRED_END_2 1

enum read_fragment_type
{
    FRAGMENT_END = 1,
    FULL_GENOME_FRAGMENT = 2,
    FULL_TRANSCRIPTOME_FRAGMENT = 3
};

enum paired_reads_strand
{
    READ_PAIRS_ARE_SAME_STRAND = 1,
    READ_PAIRS_ARE_OPP_STRAND = 2,
    READ_PAIRS_STRAND_UNKNOWN = 3
};

struct pos_in_template {
    int pos;
    int number_of_reads_in_template;
};

struct read_subtemplate {
    char* char_seq;
    char* error_str;
    
    struct penalty_array_t* fwd_penalty_array;
    struct penalty_array_t* rev_penalty_array;
    
    int length;
    struct pos_in_template pos_in_template;
    
    struct indexable_subtemplates* ists;
};

struct prior_read_information {
    int max_ref_insert_length;
    int max_fragment_length;
    enum read_fragment_type frag_type;
    enum paired_reads_strand paired_read_strand_info;
    enum assay_type_t assay;
};

struct read {
    readkey_t read_id;
    char* name;
    struct prior_read_information prior;

    struct read_subtemplate* subtemplates;
    int num_subtemplates;
};

void
init_read(
        struct read** r,
        char* readname
    );

void
free_read( struct read* r );

void
free_read_subtemplate( struct read_subtemplate* st );

void
fprintf_read( FILE* fp, struct read* r );

void
fprintf_read_subtemplate_to_fastq(
        FILE* fp,
        char* name,
        struct read_subtemplate* st
    );

/* 
   if the next readkey would be greater than maqx readkey, then
   dont return anything. negative values indicate that this should
   be ignored 
*/
int
get_next_read_from_rawread_db( 
        struct rawread_db_t* rdb,
        struct read** r,
        long max_readkey
    );

/*** Indexable subtemplates ***/
/*
 * Subsequences from read subtemplates to pass directly to the index
 */
struct indexable_subtemplate
{
    /* the offset in the read subtemplate that this index probe refers to */
    int subseq_length;
    int subseq_offset;

    /* These point into the strings that were allocated for the read
     * subtemplate */
    char* char_seq;

    /* pointers into penalty_array_t->array */
    /* these have already accounted for the subseq_offset */
    struct penalty_array_t fwd_penalty_array;
    struct penalty_array_t rev_penalty_array;
};

struct indexable_subtemplates
{
    struct indexable_subtemplate* container;
    int length;
};

void
init_indexable_subtemplate(
        struct indexable_subtemplate** ist,
        struct read_subtemplate* rst,

        int subseq_length,
        int subseq_offset
    );

void
init_indexable_subtemplates(
        struct indexable_subtemplates** ists
    );

void
free_indexable_subtemplate(
        struct indexable_subtemplate* ist
    );

void
free_indexable_subtemplates(
        struct indexable_subtemplates* ists
    );

void
add_indexable_subtemplate_to_indexable_subtemplates(
        struct indexable_subtemplate* ist,
        struct indexable_subtemplates* ists
    );

void
cache_penalty_arrays_in_read_subtemplates(
    struct read* r, struct mapping_params* mapping_params 
    );

#endif // READ_HEADER
