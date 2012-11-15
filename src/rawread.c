/* Copyright (c) 2009-2012, Nathan Boley */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <pthread.h>

#include "config.h"
#include "statmap.h"
#include "error_correction.h"
#include "rawread.h"
#include "quality.h"
#include "dna_sequence.h"
#include "util.h"
#include "read.h"

#include "log.h"

/***** Rawread initialization and utility functions *****/

void 
init_rawread( struct rawread** r,
              size_t read_len,
              size_t readname_len )
{
    *r = malloc( sizeof( struct rawread ) );

    (*r)->name = calloc( readname_len+1, sizeof(char) );
    (*r)->length = read_len;
    (*r)->char_seq = malloc( read_len*sizeof(char) );
    (*r)->error_str = malloc( read_len*sizeof(char) );
}

inline void 
free_rawread( struct rawread* r )
{
    /* free dynamically allocated strings */
    free( r->name );
    free( r->char_seq );
    free( r->error_str );

    free( r );
}

/***** FASTQ file parsing *****/

/*
 *   Move the filepointer until we see a @, indicating the start of the next 
 *   read.
 */
void
move_fastq_fp_to_next_read( FILE* fp )
{
    /* The next read begins with "\n@" */
    /* XXX - an error string could begin with @, which would also satisfy */
    long current_pos = ftell( fp );
    while( !feof(fp) ) {
        char character = fgetc( fp );
        // XXX - assumes @ only appears at the beginning of a readname
        // it can also be part of an error string, and there's no reason it
        // couldn't be used again in a read name
        if( character == '@' )
        {
            fseek( fp, current_pos, SEEK_SET );
            return;
        }
        current_pos += 1;
    }
    
    return;
}

/* 
   Reads the next line from fp, and strips trailing whitespace. 

   If 'expected_first_char' != ''. them strip the first character and 
   error if it's not expected.

   returns 0 on success, 1 on error.
*/
int safe_get_next_line( FILE* input_file, 
                        char* buffer, int buffer_size, 
                        char expected_first_char,
                        enum bool allow_empty_lines
    )
{
    char* return_code;
    
    /* get the parsed name string */
    if( '\0' != expected_first_char ) {
        char next_char = fgetc( input_file );
        if( next_char != expected_first_char )
        {
            if( feof( input_file ) )
                return 1;
            
            statmap_log( LOG_ERROR,
                    "Expected '%c' and get '%c' at position '%li'. Skipping read.",
                    expected_first_char, next_char, ftell(input_file)
                );
            move_fastq_fp_to_next_read( input_file );
            return 1;
        }
    }

    /* get the rest of the line */
    return_code = fgets( buffer, buffer_size, input_file );
    if( NULL == return_code )
    {
        if( feof(input_file) ) {
            return 1;
        } else {
            statmap_log( LOG_FATAL, "Error reading from rawread file." );
            exit( 1 );
        }
    }
        
    int buffer_len = strlen( buffer );
    
    if( buffer[buffer_len-1] != '\n' )
    {
        statmap_log( LOG_ERROR, "Unable to read full line. Skipping read." );
        move_fastq_fp_to_next_read( input_file );
        return 1;        
    }
    
    while( buffer_len > 0 && !isgraph( buffer[buffer_len-1] ) )
    {
        buffer[buffer_len-1] = '\0';
        buffer_len -= 1;
    }

    if( !allow_empty_lines && buffer_len == 0 ) {
        statmap_log( LOG_ERROR, "Line had 0 chracters. Skipping read." );
        move_fastq_fp_to_next_read( input_file );
        return 1;        
    }
    
    return 0;
}

/*
 * Populate a rawread object from a file pointer. 
 *
 * Given a file pointer (or pointers) that is (are) positioned at the begining
 * of a read, this initializes and populates the read.
 *
 * Returns 0 on success, negative values on failure. The value depends
 * upon the failure mode, read below for more info.
 *
 */
int
populate_rawread_from_fastq_file(
        FILE* input_file,
        struct rawread** r,
        enum READ_END end
    )
{
    /* Store the return of the scanf's */
    int rv;
    
    char readname[READ_BUFFER_SIZE];
    char readname2[READ_BUFFER_SIZE];
    char read[READ_BUFFER_SIZE];
    char quality[READ_BUFFER_SIZE];
    
    /***** Get the read informnation from the fastq file *****/
    while( true )
    {    
        if( feof( input_file ) )
        {
            *r = NULL;
            return EOF;
        }

        /*** Determine the read name ***/
        rv = safe_get_next_line(
            input_file, readname, READ_BUFFER_SIZE, '@', false);
        if( 1 == rv )
            continue;
        
        /*
         * Since we aren't relying on the /1 and /2 suffixes to determine
         * read end anymore, end should always be known
         */
        assert( end != UNKNOWN );

        /* If the read end is known, eliminate the slash */
        if( end == FIRST || end == SECOND )
        {
            /* find the final / */
            char* fwd_slash_pntr = strrchr ( readname, '/' );
            /* set the slash to '\0', to eliminate it from the read name */
            if( NULL != fwd_slash_pntr )
                *fwd_slash_pntr = '\0';
        } else if( end == NORMAL )   {
            /* if the end is normal, we dont need to do antyhing */
        } else {
            statmap_log( LOG_ERROR, "Read end should be know" );
            assert(false);
        }
            
        /*** get the actual read */
        rv = safe_get_next_line(input_file, read, READ_BUFFER_SIZE, '\0', false);
        if( 1 == rv )
            continue;
        
        /*** get the next read name - we discard this */
        rv = safe_get_next_line(
            input_file, readname2, READ_BUFFER_SIZE, '+', true);
        if( 1 == rv )
            continue;
        
        /*** get the quality score */
        rv = safe_get_next_line(
            input_file, quality, READ_BUFFER_SIZE, '\0', false);
        if( 1 == rv )
            continue;
        
        /* if everything has been loaded correctly, then break out 
           of the while loop, and store the info in the structure */
        break;
    }

    /* calculate the length of the read and the readname */
    int readname_len = strlen( readname );
    int read_len = strlen( read );
    
    /***** Initialize and Populate the raw read *****/ 
    init_rawread( r, read_len, readname_len );

    /* Already set by the init function */
    assert( (*r)->length == read_len );
    
    /* Copy the readname */
    memcpy( (*r)->name, readname, sizeof(char)*(readname_len + 1) );
    
    /* Copy the read */
    memcpy( (*r)->char_seq, read, sizeof(char)*(read_len) );
    
    /* Copy the error string */
    memcpy( (*r)->error_str, quality, sizeof(char)*(read_len) );    

    /* Set the read end */
    (*r)->end = end;
    
    return 0;
}

/******* BEGIN raw read DB code ********************************************/

void 
init_rawread_db( struct rawread_db_t** rdb )
{
    (*rdb) = malloc( sizeof( struct rawread_db_t ) );

    (*rdb)->readkey = 0;

    pthread_mutexattr_t mta;
    pthread_mutexattr_init(&mta);
    (*rdb)->mutex = malloc( sizeof(pthread_mutex_t) );
    pthread_mutex_init( (*rdb)->mutex, &mta );

    (*rdb)->single_end_reads = NULL;
    (*rdb)->paired_end_1_reads = NULL;
    (*rdb)->paired_end_2_reads = NULL;

    (*rdb)->unmappable_single_end_reads = NULL;
    (*rdb)->unmappable_paired_end_1_reads = NULL;
    (*rdb)->unmappable_paired_end_2_reads = NULL;

    (*rdb)->non_mapping_single_end_reads = NULL;
    (*rdb)->non_mapping_paired_end_1_reads = NULL;
    (*rdb)->non_mapping_paired_end_2_reads = NULL;

    (*rdb)->file_type = UNKNOWN;
    return;
}

void 
close_rawread_db( struct rawread_db_t* rdb )
{
    if( rdb->single_end_reads != NULL )
    {
        fclose( rdb->single_end_reads );
        fclose( rdb->unmappable_single_end_reads );
        fclose( rdb->non_mapping_single_end_reads );
    } else  {
        fclose( rdb->paired_end_1_reads );
        fclose( rdb->unmappable_paired_end_1_reads );
        fclose( rdb->non_mapping_paired_end_1_reads );

        fclose( rdb->paired_end_2_reads );
        fclose( rdb->unmappable_paired_end_2_reads );
        fclose( rdb->non_mapping_paired_end_2_reads );
    }
    
    pthread_mutex_destroy( rdb->mutex );
    free( rdb->mutex );
    
    free( rdb );
    return;
}

void
lock_rawread_db(
        struct rawread_db_t* rdb
    )
{
    //fprintf(stderr, "BEFORE lock on rk %d, PID : %u\n", rdb->readkey, pthread_self());
    pthread_mutex_lock( rdb->mutex );
    //fprintf(stderr, "Locked on readkey %d, PID : %u\n", rdb->readkey, pthread_self() );
}

void
unlock_rawread_db(
        struct rawread_db_t* rdb
    )
{
    //fprintf(stderr, "Unlocked on readkey %d, PID : %u\n", rdb->readkey, pthread_self() );
    pthread_mutex_unlock( rdb->mutex );
}

void
add_single_end_reads_to_rawread_db( 
    struct rawread_db_t* rdb, char* rifname, 
    enum inputfile_type iftype, enum assay_type_t assay )
{
    /* Make sure that no files have been added yet */
    assert( rdb->paired_end_1_reads == NULL );
    assert( rdb->paired_end_2_reads == NULL );
    assert( rdb->single_end_reads == NULL );

    char buffer[500];
    assert( strlen( rifname ) < 450 );

    /* Pair 1 */
    rdb->single_end_reads = open_check_error( rifname, "r" );

    strcpy( buffer, rifname );
    strcat( buffer, ".unmappable" );
    rdb->unmappable_single_end_reads = open_check_error( buffer, "a" );
    
    strcpy( buffer, rifname );
    strcat( buffer, ".nonmapping" );
    rdb->non_mapping_single_end_reads = open_check_error( buffer, "a" );
    
    /* For now, we only support fastq files */
    assert( iftype == FASTQ );
    rdb->file_type = iftype;
    
    rdb->assay = assay;
    
    return;
}

void
add_paired_end_reads_to_rawread_db( 
    struct rawread_db_t* rdb,
    char* rifname1, char* rifname2,
    enum inputfile_type iftype,
    enum assay_type_t assay )
{
    char buffer[500];
    assert( strlen( rifname1 ) < 480 );
    assert( strlen( rifname2 ) < 480 );

    /* Make sure that no files have been added yet */
    assert( rdb->paired_end_1_reads == NULL );
    assert( rdb->paired_end_2_reads == NULL );
    assert( rdb->single_end_reads == NULL );
    
    /* Pair 1 */
    rdb->paired_end_1_reads = open_check_error( rifname1, "r" );
    
    strcpy( buffer, rifname1 );
    strcat( buffer, ".unmappable" );
    rdb->unmappable_paired_end_1_reads = open_check_error( buffer, "w" );
    
    strcpy( buffer, rifname1 );
    strcat( buffer, ".nonmapping" );
    rdb->non_mapping_paired_end_1_reads = open_check_error( buffer, "w" );
    
    /* Pair 2 */
    rdb->paired_end_2_reads = open_check_error( rifname2, "r" );
    
    strcpy( buffer, rifname2 );
    strcat( buffer, ".unmappable" );
    rdb->unmappable_paired_end_2_reads = open_check_error( buffer, "w" );
    
    strcpy( buffer, rifname2 );
    strcat( buffer, ".nonmapping" );
    rdb->non_mapping_paired_end_2_reads = open_check_error( buffer, "w" );
    
    /* File Type */
    assert( iftype == FASTQ );
    rdb->file_type = iftype;
    
    rdb->assay = assay;
    
    return;    
}

void
rewind_rawread_db( struct rawread_db_t* rdb )
{
    /* reset the read key */
    rdb->readkey = 0;

    /* rewind all of the file pointers */

    if( rdb->paired_end_1_reads != NULL )
        rewind( rdb->paired_end_1_reads );
    
    if( rdb->paired_end_2_reads != NULL )
        rewind( rdb->paired_end_2_reads );
    
    if( rdb->single_end_reads != NULL )
        rewind( rdb->single_end_reads );
    
    return;
}

enum bool
rawread_db_is_empty( struct rawread_db_t* rdb )
{
    if( rdb->single_end_reads != NULL )
    {
        if( feof( rdb->single_end_reads ) )
            return true;
        return false;
    } 
    /* If the reads are paired */
    else {
        assert( rdb->paired_end_1_reads != NULL );
        assert( rdb->paired_end_2_reads != NULL );
        if( feof( rdb->paired_end_1_reads ) )
        {
            if( !feof(rdb->paired_end_2_reads) ) 
            {
                statmap_log( LOG_ERROR, "File Mismatch: Paired files have different numbers of reads" );
                exit( -1 );
            }
            return true;
        }
        return false;
    }        
    
}

struct rawread_db_state
save_rawread_db_state(
        struct rawread_db_t *rdb
    )
{
    assert( rdb->mutex != NULL );
    pthread_mutex_lock( rdb->mutex );

    /* Since the rawread db stores *either* single or paired end reads, it
     * should be the case that either the single_end_reads fp has been
     * initialized, or both paired end read fps have been initialized */
    assert( rdb->single_end_reads != NULL || (
                rdb->single_end_reads == NULL && (
                    rdb->paired_end_1_reads != NULL &&
                    rdb->paired_end_2_reads != NULL )));

    struct rawread_db_state saved_state;
    saved_state.readkey = rdb->readkey;

    /* Initialize the read positions to -1 for "unset" */
    saved_state.single_end_reads_pos = -1;
    saved_state.paired_end_1_reads_pos = -1;
    saved_state.paired_end_2_reads_pos = -1;

    if( rdb->single_end_reads != NULL ) {
        saved_state.single_end_reads_pos = ftell(rdb->single_end_reads);
        assert( saved_state.single_end_reads_pos != -1 );
    } else {
        saved_state.paired_end_1_reads_pos = ftell(rdb->paired_end_1_reads);
        assert( saved_state.paired_end_1_reads_pos != -1 );
        saved_state.paired_end_2_reads_pos = ftell(rdb->paired_end_2_reads);
        assert( saved_state.paired_end_2_reads_pos != -1 );
    }

    pthread_mutex_unlock( rdb->mutex );

    return saved_state;
}

void
restore_rawread_db_state(
        struct rawread_db_t* rdb,
        struct rawread_db_state state
    )
{
    assert( rdb->mutex != NULL );
    pthread_mutex_lock( rdb->mutex );

    /* reset the readkey */
    rdb->readkey = state.readkey;

    int rv;
    /* reset the underlying file pointers */
    if( state.single_end_reads_pos != -1 ) {
        rv = fseek( rdb->single_end_reads, state.single_end_reads_pos,
                SEEK_SET );
        assert( rv != -1 );
    } else {
        rv = fseek( rdb->paired_end_1_reads, state.paired_end_1_reads_pos,
                SEEK_SET );
        assert( rv != -1 );
        rv = fseek( rdb->paired_end_2_reads, state.paired_end_2_reads_pos,
                SEEK_SET );
        assert( rv != -1 );
    }

    pthread_mutex_unlock( rdb->mutex );
}

/******* END raw read DB code ********************************************/
