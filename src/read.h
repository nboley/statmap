/* Copyright (c) 2009-2010 Nathan Boley */

#ifndef READS_HEADER
#define READS_HEADER

#include "config.h"

struct read_subtemplate {
    char* char_seq;
    char* error_seq;
    int length;
    int pos_in_template;    // +/- indexed, like Python list indexing
    enum READ_END end;
};

struct read {
    char* name;
    enum assay_type_t assay;

    struct read_subtemplate* subtemplates;
    int num_subtemplates;
};

void
init_read( struct read** r,
           size_t readname_len );

void
free_read( struct read* r );

void
add_subtemplate_to_read( struct read* r,
        char* char_seq, char* error_seq,
        int len, int pos_in_template,
        enum READ_END end
    );

void
fprintf_read( FILE* fp, struct read* r );

/* determine whether reads are mappable */
enum bool
filter_read( struct read* r,
             struct error_model_t* error_model );

/* 
   if the next readkey would be greater than maqx readkey, then
   dont return anything. negative values indicate that this should
   be ignored 
*/
int
get_next_read_from_rawread_db( 
    struct rawread_db_t* rdb, readkey_t* readkey,
    struct read* r,
    long max_readkey );

#endif // READS_HEADER
