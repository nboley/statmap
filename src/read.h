/* Copyright (c) 2009-2012 Nathan Boley */

#ifndef READ_HEADER
#define READ_HEADER

#include "config.h"
#include "rawread.h"
#include "quality.h"

//struct penalty_array_t; // fwd declaration (?)

#define POS_SINGLE_END 1
#define POS_PAIRED_END_1 1
#define POS_PAIRED_END_2 -1

struct pos_in_template {
    int pos;
    int number_of_reads_in_template;
    enum bool is_full_fragment;
};

struct read_subtemplate {
    char* char_seq;
    char* error_str;
    int length;
    struct pos_in_template pos_in_template;
};

struct read {
    char* name;
    enum assay_type_t assay;

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

/* determine whether reads are mappable */
enum bool
filter_read(
        struct read* r,
        struct error_model_t* error_model
    );

/* 
   if the next readkey would be greater than maqx readkey, then
   dont return anything. negative values indicate that this should
   be ignored 
*/
int
get_next_read_from_rawread_db( 
        struct rawread_db_t* rdb,
        readkey_t* readkey,
        struct read** r,
        long max_readkey
    );

/*** Indexable subtemplates ***/
/*
 * Subsequences from read subtemplates to pass directly to the index
 */
struct indexable_subtemplate
{
    int subseq_offset;

    /* These point into the strings that were allocated for the read
     * subtemplate */
    char* char_seq;

    /* pointers into penalty_array_t->array */
    struct penalty_t* fwd_penalties;
    struct penalty_t* rev_penalties;
};

struct indexable_subtemplates
{
    struct indexable_subtemplate* container;
    int length;
};

void
init_indexable_subtemplate(
        struct indexable_subtemplate** ist,

        int subseq_offset,
        char* char_seq,

        struct penalty_array_t* fwd_penalty_array,
        struct penalty_array_t* rev_penalty_array
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
