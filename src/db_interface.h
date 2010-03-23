/* Copyright (c) 2009-2010 Nathan Boley */

#ifndef DB_INTERFACE
#define DB_INTERFACE

#include "mapped_location.h"

/* 
 * store the data structures necessary to write a read 
 *
 * this is a generic interface so that we can store reads in
 * multiple locations - including flat files, and a berkeley DB
 * database for now but, eventually, we'd like to cupport postgres.
 *
 * Add/Get candidate_mappings through interface - namely
 *
 * add_candidate_mappings_to_db
 * get_candidate_mapping_from_db
 *
 */
typedef struct {
    /* 
     * pointer to a berkeley db database 
     * NULL if unused
     */
    void* bdb; 
    
    /* eventually, a pointer to a postgres DB */
    // PGDB* pdb;

    /* The directory that stores all of the data */
    char* data_dir;
    
    /* buffer the reads for insert */
    candidate_mappings* buffer;

    /* pointers to flat files that store candidate mappings */
    /* NULL indicates that the fp is not used */
    /* full reads - ie non partial reads */

    /***** Files to store 'normal' reads ********************/
    /** store non-paired end reads **/
    /*** the fwd strand reads ***/
    FILE** full_fwd_rds;
    /*** the bkwd strand reads ***/
    FILE** full_bkwd_rds;

    /** store paired end reads **/
    /*** the first paired end  ***/
    /**** the fwd strand reads ****/
    FILE** full_1st_fwd_rds;
    /**** the bkwd strand reads ****/
    FILE** full_1st_bkwd_rds;
    /*** the second paired end  ***/
    /**** the fwd strand reads ****/
    FILE** full_2nd_fwd_rds;
    /**** the bkwd strand reads ****/
    FILE** full_2nd_bkwd_rds;

    /***** Files to store 'junction' reads ********************/
    /* BIG TODO */

} candidate_mappings_db;

int
cmp_joining_queue_datum( const void* a, const void* b );

#define MAX_KEY_SIZE 254

typedef struct {
    long read_id;
    candidate_mapping mapping;
    FILE* stream; 
} joining_queue_datum;

typedef struct {
    candidate_mappings_db* db;
    joining_queue_datum** mappings_queue;
    int queue_len;
    long curr_read_id;
} candidate_mappings_db_cursor;



void
init_candidate_mappings_db( candidate_mappings_db* db, 
                            char* fname_prefix );

void
close_candidate_mappings_db( candidate_mappings_db* db );

void
add_candidate_mappings_to_db( 
    candidate_mappings_db* db, 
    candidate_mappings* mappings,
    long read_id,
    int thread_num );

#define CURSOR_EMPTY 1

void
open_candidate_mappings_cursor(
    candidate_mappings_db* db,
    candidate_mappings_db_cursor** cursor);

void
close_candidate_mappings_cursor( 
    candidate_mappings_db_cursor* cursor );


int
get_next_candidate_mapping_from_cursor( 
    candidate_mappings_db_cursor* cursor, 
    candidate_mappings** mappings,
    long* read_id );

void
join_all_candidate_mappings( 
    candidate_mappings_db* cand_mappings_db,
    mapped_reads_db* mpd_rds_db );


#endif // DB_INTERFACE define


