/* Copyright (c) 2009-2012 Nathan Boley */

#ifndef READ_HEADER
#define READ_HEADER

#include "config.h"
#include "rawread.h"
#include "quality.h"

struct penalty_array_t; // fwd declaration (?)

#define POS_SINGLE_END 0
#define POS_PAIRED_END_1 0
#define POS_PAIRED_END_2 1

struct pos_in_template {
    int pos;
    int number_of_reads_in_template;
    enum bool is_full_fragment;
};

struct read_subtemplate {
    char* char_seq;
    char* error_str;
    
    struct penalty_array_t* fwd_penalty_array;
    struct penalty_array_t* rev_penalty_array;
    
    int length;
    struct pos_in_template pos_in_template;
};

struct prior_read_information {
    int max_ref_insert_length;
    int max_fragment_length;
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
    struct penalty_t* fwd_penalties;
    struct penalty_t* rev_penalties;

    /* expected value of this index probe */
    float expected_value;
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

#endif // READ_HEADER
