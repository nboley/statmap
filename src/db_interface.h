/* Copyright (c) 2009-2010 Nathan Boley */

#ifndef DB_INTERFACE
#define DB_INTERFACE

#include "mapped_location.h"
// #include "fastq.h"

/* 
 * Usages of the database interface
 *
 *
 * Temporary Mapped Read Storage
 *
 * 1) Read a read in ( currently, with fread via get_next_read )
 *    i.   Mark the strand, paired end side, ( that's all we know )
 *    ii.  Map it
 *         a. For each mapping create a result
 *            NOTE: The result is of type mapped_location 
 *                  which can be found in mapped_location.h
 *            1. Mark the chromosome
 *            2. Mark the mapped strand
 *            3. Mark the side of the paired end
 *            4. Mark the bp position
 *            5. Mark the penalty ( as a log probability )
 *            6. If this is mrna seq
 *               i.   Mark the gene direction
 *                    a. for normal reads, add both gene directions
 *               ii.  Mark the read type
 *                    a. 'Normal' Read
 *                    b.  Junction Read
 *                        NOTE: Currently, this is part of genome location type, 
 *                              but in the future this may be a recheck
 *                        1. Mark as a junction read
 *                        2. Mark pre-intron, post-intron, or novel ( maybe )
 *                    c.  Polya Read
 *         b. Add the mapped result to the result set
 *            1. Use the add_read ( or add_reads_from_mapped_locations ) in 
 *               the db_interface
 *               NOTE: We consider a read to be the entire paired end. That is
 *                     for UNIQUENAME/1 is just the read key UN with read 1. 
 *               TODO: Make the sort order of these better, such that they use
 *                     chr/gene_strnd/read_type/paired_end_side/penalty
 *    iii. Merge the mapped result set for each mapped read type
 *         TODO: How do we know when every read of this type has been mapped? 
 *         A separate db? For paired end reads, there are 4. For single, 2. ( One
 *         for each potential read direction ) 
 *         a. Group possible reads
 *            - group by chr
 *                - group by gene strnd
 *                    - group by read type 
 *         b. Pair junction reads
 *
 */

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
    char key[MAX_KEY_SIZE+1];
    unsigned long read_id;
    candidate_mapping mapping;
    FILE* stream; 
} joining_queue_datum;

typedef struct {
    candidate_mappings_db* db;
    joining_queue_datum** mappings_queue;
    int queue_len;
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
    int thread_num,
    char* readname,
    unsigned long read_id);

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
    char* key );



#endif // DB_INTERFACE define


