/* Copyright (c) 2009-2012 Nathan Boley */

#ifndef RAWREADS_HEADER
#define RAWREADS_HEADER

#include <math.h>
#include <pthread.h>

#include "config.h"

#define READ_BUFFER_SIZE 65536

typedef unsigned int readkey_t;

/* 
 * Define a struct to hold read data. 
 * 
 * Every read that we get fom a platform specific source should be 
 * converted into one of these, and then written to the pre-processed db.
 * Then, we map from the precprocessed db.
 *
 */

/*
 * Exact representation of read information stored in FASTA file
 */
struct rawread {
    char* name;
    int length;
    char* char_seq;
    char* error_str;
    enum READ_END end;
};

/* fwd declaration for error_model_t */
struct error_model_t;

/* Initialize a raw read. The length's are needed to init the char strings */
void 
init_rawread( struct rawread** r,
              size_t read_len,
              size_t readname_len );

/* Free a raw read. */
void 
free_rawread( struct rawread* r );

/* Populate a read from the next read in a fastq file */

int
populate_rawread_from_fastq_file(
        FILE* input_file,
        struct rawread** r,
        enum READ_END end
    );

/**************** Raw Read DB **********************/
/* An API for consolidating the many types of 
 * reads that we may encounter. This deals with
 * read order, thread locking, buffered read 
 * optimizations, etc. 
 *
 * The basic approach for usage is to
 * 1) Initialize a DB ( init_rawread_db )
 * 2) Add flat files to the DB ( add_reads_to_rawread_db )
 * 
 * When all of the reads are added
 *
 * 1) Call get_next_read_from_rawread_db() 
 *     until rawread_db_is_empty() === TRUE.
 * 
 * 2) rewind_rawread_db() to move the pointer to the begginging
 * 3) call get_next_read to join all of the reads
 *
 */

enum inputfile_type {
    /* UNKNOWN = 0 */
    FASTQ = 1
};

struct rawread_db_t {
    /* 
       This 'db' is functioning as both a cursor and a db. Therefore, 
       readkey stores the key of the current read. Calling get_next_read
       will pull out the next read, and increment this. rewind will
       reset this to 0. 
    */
    readkey_t readkey;
    /* We lock this whenever we grab a read */
    //pthread_spinlock_t* lock; 
    pthread_mutex_t* mutex;

    FILE* single_end_reads;
    FILE* paired_end_1_reads;
    FILE* paired_end_2_reads;

    FILE* unmappable_single_end_reads;
    FILE* unmappable_paired_end_1_reads;
    FILE* unmappable_paired_end_2_reads;

    FILE* non_mapping_single_end_reads;
    FILE* non_mapping_paired_end_1_reads;
    FILE* non_mapping_paired_end_2_reads;

    enum inputfile_type file_type;
    enum assay_type_t assay;
};

void 
init_rawread_db( struct rawread_db_t** rdb );

void 
close_rawread_db( struct rawread_db_t* rdb );

void
lock_rawread_db(
        struct rawread_db_t* rdb
    );

void
unlock_rawread_db(
        struct rawread_db_t* rdb
    );

void
add_single_end_reads_to_rawread_db( 
    struct rawread_db_t* rdb, char* rifname, 
    enum inputfile_type iftype,
    enum assay_type_t assay );

void
add_paired_end_reads_to_rawread_db( 
    struct rawread_db_t* rdb,
    char* rifname1, char* rifname2,
    enum inputfile_type iftype,
    enum assay_type_t assay );

void
rewind_rawread_db( struct rawread_db_t* rdb );

enum bool
rawread_db_is_empty( struct rawread_db_t* rdb );

void
move_fp_to_next_read( FILE* fp );

/**************** END Raw Read DB **********************/

#endif // RAWREADS_HEADER
