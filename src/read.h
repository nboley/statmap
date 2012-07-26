/* Copyright (c) 2009-2010 Nathan Boley */

#ifndef READ_HEADER
#define READ_HEADER

#include "config.h"
#include "rawread.h"

struct subtemplate {
    char* char_seq;
    char* error_str;
    int length;
    int pos_in_template;    // +/- indexed, like Python list indexing
    enum READ_END end;
};

struct read {
    char* name;
    enum assay_type_t assay;

    struct subtemplate* r1;
    struct subtemplate* r2;
};

void
init_read(
        struct read** r,
        char* readname
    );

void
free_read( struct read* r );

void
free_subtemplate( struct subtemplate* st );

void
fprintf_read( FILE* fp, struct read* r );

void
fprintf_subtemplate_to_fastq( FILE* fp, char* name, struct subtemplate* st );

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

#endif // READ_HEADER
