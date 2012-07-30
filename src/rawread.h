/* Copyright (c) 2009-2012 Nathan Boley */

#ifndef READS_HEADER
#define READS_HEADER

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

struct rawread {
    /* should equal to strlen( char_seq ) */
    unsigned char length;
    char* name;
    char* char_seq;
    char* error_str;
    enum READ_END end;
    enum STRAND strand;
    enum assay_type_t assay;
    // LETTER_TYPE* packed_seq;
    // float* penalties;
    // float* inverse_penalties; 
    // float* mutation_type_penalties;
};

/* fwd declaration for error_model_t */
struct error_model_t;

/* Initialize a raw read. The length's are needed to init the char strings */
void 
init_rawread( struct rawread** r,
              int seq_len,
              size_t readname_len);

/* Free a raw read. */
void 
free_rawread( struct rawread* r );

void
fprintf_rawread( FILE* fp, struct rawread* r );

void
fprintf_rawread_to_fastq( FILE* fastq_fp, struct rawread* r );

void
marshal_rawread( struct rawread* r, char** buffer, size_t* buffer_size );

void
unmarshal_rawread( struct rawread** r, char* buffer );

/* Populate a read from the next read in a fastq file */

int
populate_read_from_fastq_file( FILE* f, struct rawread** r, enum READ_END end );

/* determine whether reads are mappable */
enum bool
filter_rawread( struct rawread* r,
                struct error_model_t* error_model );

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
    /* We lock this mutex whenever we grab a read */
    pthread_spinlock_t* lock; 

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

/* 
   if the next readkey would be greater than maqx readkey, then
   dont return anything. negative values indicate that this should
   be ignored 
*/
int
get_next_read_from_rawread_db( 
    struct rawread_db_t* rdb, readkey_t* readkey,
    struct rawread** r1, struct rawread** r2,
    long max_readkey );

/**************** END Raw Read DB **********************/

#endif // READS_HEADER


