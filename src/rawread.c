/* Copyright (c) 2009-2010, Nathan Boley */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <pthread.h>
#include <assert.h>

#include "config.h"
#include "statmap.h"
#include "error_correction.h"
#include "rawread.h"
#include "quality.h"
#include "dna_sequence.h"

static FILE* 
open_check_error( char* fname, char* file_mode )
{
    FILE* tmp;
    tmp = fopen( fname, file_mode );
    if( tmp == NULL )
    {
        fprintf( stderr, "Error opening '%s\n'", fname );
        exit( -1 );
    }
    return tmp;
}

inline void 
init_rawread( struct rawread** r,
              int seq_len,
              size_t readname_len  )
{
    *r = malloc( sizeof( struct rawread ) );

    /* BUG TODO check for nulls */
    (*r)->length = seq_len;
    (*r)->name = malloc( sizeof(char)*(readname_len+1) );
    (*r)->char_seq = malloc( sizeof(char)*seq_len );
    (*r)->error_str = malloc( sizeof(char)*seq_len );
}

inline void 
free_rawread( struct rawread* r )
{
    if( r->name != NULL )
        free( r->name );

    if( r->char_seq != NULL )    
        free( r->char_seq );

    if( r->error_str != NULL )
        free( r->error_str );
    
    free( r );
}

void
fprintf_rawread( FILE* fp, struct rawread* r )
{
    fprintf(fp, "%u\n", r->length);
    fprintf(fp, "%s\n", r->name);
    fprintf(fp, "%.*s\n", r->length, r->char_seq);
    fprintf(fp, "%.*s\n", r->length, r->error_str);
    fprintf(fp, "%u\n", r->end);
    fprintf(fp, "%u\n", r->strand);
    fprintf(fp, "\n\n");
}

void
fprintf_rawread_to_fastq( FILE* fastq_fp, struct rawread* r )
{
    
    fprintf(fastq_fp, "@%s\n", r->name);
    fprintf(fastq_fp, "%.*s\n", r->length, r->char_seq);    
    fprintf(fastq_fp, "+%s\n", r->name);
    fprintf(fastq_fp, "%.*s\n", r->length, r->error_str);
}
/*
 *   Move the filepointer until we see a @, indicating the start of the next 
 *   read.
 */
void
move_fastq_fp_to_next_read( FILE* fp )
{
    long current_pos = ftell( fp );
    while( !feof(fp) ) {
        char character = fgetc( fp );
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
 * Populate a rawread object from a file pointer. 
 *
 * Given a file pointer that is positioned at the begining of a read
 * this initializes and populates a read.
 *
 * Returns 0 on success, negative values on failure. The value depends
 * upon the failure mode, read below for more info.
 *
 */

enum READ_END
determine_read_end_from_readname( char* readname )
{
    /* find the final / */
    char* fwd_slash_pntr = 
        strrchr ( readname, '/' );
    
    enum READ_END end = UNKNOWN;
    
    /* if there is no slash, assume its not a paired end read */
    if( fwd_slash_pntr == NULL )
    {
        return NORMAL;
    } else {
        /* Determine the read end from integer following the slash */
        switch( atoi( fwd_slash_pntr + 1 ) )
        {
        case 1:
            end = FIRST;
            break;
        case 2:
            end = SECOND;
            break;
        default:
            fprintf(stderr, 
                    "Unrecognized read end '%i'\n", 
                    atoi( fwd_slash_pntr + 1 ) );
            exit( -1 );
        }

        /* set the slash to '\0', to eliminate it from the read name */
        if( NULL != fwd_slash_pntr )
            *fwd_slash_pntr = '\0';
        
        return end;
    }
    
    /* we should never get here */
    assert( false );
}

/* 
   Reads the next line from fp, and striped trailing whitespace. 

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
            
            fprintf(stderr, "ERROR:    Expected '%c' and get '%c' at position '%li'. Skipping read.\n",
                    expected_first_char, next_char, ftell(input_file) );
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
            perror( "FATAL:    Error reading from rawread file." );
            exit( 1 );
        }
    }
        
    int buffer_len = strlen( buffer );
    
    if( buffer[buffer_len-1] != '\n' )
    {
        fprintf(stderr, "ERROR:    Unable to read full line. Skipping read.\n");
        move_fastq_fp_to_next_read( input_file );
        return 1;        
    }
    
    while( buffer_len > 1 && !isgraph( buffer[buffer_len-1] ) )
    {
        buffer[buffer_len-1] = '\0';
        buffer_len -= 1;
    }

    if( !allow_empty_lines && buffer_len == 0 ) {
        fprintf(stderr, "ERROR:    Line had 0 chracters. Skipping read.\n");
        move_fastq_fp_to_next_read( input_file );
        return 1;        
    }
    
    return 0;
}

int
populate_read_from_fastq_file( 
    FILE* input_file, struct rawread** r, enum READ_END end )
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
            assert( end == UNKNOWN );
            end = determine_read_end_from_readname( readname );
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
    
    /* Set the strand ( FIXME ) when this is known */
    (*r)->strand = UNKNOWN;

    return 0;
}

enum bool
filter_rawread( struct rawread* r,
                struct error_model_t* error_model )
{
    /* Might pass a NULL read (r2 in the single read case) 
       from find_candidate_mappings */
    if( r == NULL )
        return false; // Do not filter nonexistent read

    /***************************************************************
     * check to make sure this read is mappable
     * we consider a read 'mappable' if:
     * 1) There are enough hq bps
     *
     */

    /* Make sure the global option has been set 
       ( it's initialized to -1 ); */
    assert( min_num_hq_bps >= 0 );

    int num_hq_bps = 0;
    int i;
    for( i = 0; i < r->length; i++ )
    {
        /*
           compute the inverse probability of error (quality)
           NOTE when error_prb receieves identical bp's, it returns the
           inverse automatically
         */
        double error = error_prb( r->char_seq[i], r->char_seq[i], 
                                  r->error_str[i], i, error_model );
        double qual = pow(10, error );
        
        /* count the number of hq basepairs */
        if( qual > 0.999 )
            num_hq_bps += 1;
    }

    if ( num_hq_bps < min_num_hq_bps )
        return true;
        
    return false;
}

/******* BEGIN raw read DB code ********************************************/

void 
init_rawread_db( struct rawread_db_t** rdb )
{
    (*rdb) = malloc( sizeof( struct rawread_db_t ) );

    (*rdb)->readkey = 0;

    (*rdb)->lock = malloc( sizeof(pthread_spinlock_t) );
    pthread_spin_init( (*rdb)->lock, PTHREAD_PROCESS_PRIVATE );

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
    
    pthread_spin_destroy( rdb->lock );
    free( (void*) rdb->lock );
    
    free( rdb );
    return;
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
                fprintf( stderr, "ERROR       :  File Mismatch: Paired files have different numbers of reads\n" );
                exit( -1 );
            }
            return true;
        }
        return false;
    }        
    
}

int
get_next_read_from_rawread_db( 
    struct rawread_db_t* rdb, readkey_t* readkey, 
    struct rawread** r1, struct rawread** r2,
    long max_readkey )
{
    pthread_spin_lock( rdb->lock );
    
    if( max_readkey >= 0
        && rdb->readkey >= (readkey_t) max_readkey )
    {
        pthread_spin_unlock( rdb->lock );
        *r1 = NULL;
        *r2 = NULL;
        return EOF;
    }

    /* 
     * Store the return value - 
     * 0 indicates success, negative failure , EOF no more reads.
     */
    int rv = -10;

    /* If the reads are single ended */
    if( rdb->single_end_reads != NULL )
    {
        /* indicate that the reads are single ended */
        *r2 = NULL;
        
        rv = populate_read_from_fastq_file( 
            rdb->single_end_reads, r1, NORMAL );
        if( rv == EOF )
        {
            *r1 = NULL;
            pthread_spin_unlock( rdb->lock );
            return EOF;
        }        

        (*r1)->assay = rdb->assay;
    } 
    /* If the reads are paired */
    else {
        assert( rdb->paired_end_1_reads != NULL );
        assert( rdb->paired_end_2_reads != NULL );
        
        rv = populate_read_from_fastq_file( 
            rdb->paired_end_1_reads, r1, FIRST );
        if ( rv == EOF )
        {
            /* Make sure the mate is empty as well */
            rv = populate_read_from_fastq_file( 
                rdb->paired_end_2_reads, r2, SECOND );
            assert( EOF == rv );
            assert( rawread_db_is_empty( rdb ) );
            
            *r1 = NULL;
            *r2 = NULL;
            pthread_spin_unlock( rdb->lock );
            return EOF;
        }
        
        rv = populate_read_from_fastq_file( 
            rdb->paired_end_2_reads, r2, SECOND );
        
        /* if returned an EOF, then it should have been returned for r1 */
        assert( rv == 0 );
        
        (*r1)->assay = rdb->assay;
        (*r2)->assay = rdb->assay;
    }

    /* increment the read counter */
    *readkey = rdb->readkey;
    rdb->readkey += 1;

    pthread_spin_unlock( rdb->lock );
    return 0;
}


/******* END raw read DB code ********************************************/
