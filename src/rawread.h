/* Copyright (c) 2009-2010 Nathan Boley */

#ifndef READS_HEADER
#define READS_HEADER

#include "config.h"

/* 
 * Define a struct to hold read data. 
 * 
 * Every read that we get fom a platform specific source should be 
 * converted into one of these, and then written to the pre-processed db.
 * Then, we map from the precprocessed db.
 *
 */

typedef struct {
    /* should equal to strlen( char_seq ) */
    unsigned char length;
    char* name;
    char* char_seq;
    char* error_str;
    enum READ_END end;
    enum STRAND strand;
    // LETTER_TYPE* packed_seq;
    // float* penalties;
    // float* inverse_penalties; 
    // float* mutation_type_penalties;
} rawread;

/* Initialize a raw read. The length's are needed to init the char strings */
void 
init_rawread( rawread** r,
              int seq_len,
              size_t readname_len);

/* Free a raw read. */
void 
free_rawread( rawread* r );

void
fprintf_rawread( FILE* fp, rawread* r );

void
fprintf_rawread_to_fastq( FILE* fastq_fp, rawread* r );

void
marshal_rawread( rawread* r, char** buffer, size_t* buffer_size );

void
unmarshal_rawread( rawread** r, char* buffer );

/* Populate a read from the next read in a fastq file */

int
populate_read_from_fastq_file( FILE* f, rawread** r );

/* determine whether reads are mappable */
enum bool
filter_rawread( rawread* r );


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

typedef struct {
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
} rawread_db_t;

void 
init_rawread_db( rawread_db_t** rdb );

void 
close_rawread_db( rawread_db_t* rdb );

void
add_single_end_reads_to_rawread_db( 
    rawread_db_t* rdb, char* rifname, 
    enum inputfile_type iftype );

void
add_paired_end_reads_to_rawread_db( 
    rawread_db_t* rdb,
    char* rifname1, char* rifname2,
    enum inputfile_type iftype );

void
rewind_rawread_db( rawread_db_t* rdb );

enum bool
rawread_db_is_empty( rawread_db_t* rdb );

int
get_next_mappable_read_from_rawread_db( 
    rawread_db_t* rdb, rawread** r1, rawread** r2 );

int
get_next_read_from_rawread_db( 
    rawread_db_t* rdb, rawread** r1, rawread** r2 );

/**************** END Raw Read DB **********************/

#endif // READS_HEADER


