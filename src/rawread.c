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
 * Populate a rawread object from a file pointer. 
 *
 * Given a file pointer that is positioned at the begining of a read
 * this initializes and populates a read.
 *
 * Returns 0 on success, negative values on failure. The value depends
 * upon the failure mode, read below for more info.
 *
 */
int
populate_read_from_fastq_file( 
    FILE* input_file, struct rawread** r, enum READ_END end )
{
    /* Store the return of the scanf's */
    int return_code;
    
    char readname[200];
    char read[200];
    char quality[200];
    
    /***** Get the read informnation from the fastq file *****/

    /*** Determine the read name ***/
    /* get the parsed name string */
    return_code = fscanf( input_file, "@%s\n", readname );

    if( return_code == EOF )
    {
        assert( feof( input_file ) );
        *r = NULL;
        return EOF;
    }
    
    /* If the read end is known, eliminate the slash */
    if( end == FIRST || end == SECOND )
    {
        /* find the final / */
        char* fwd_slash_pntr = 
            strrchr ( readname, '/' );

       /* set the slash to '\0', to eliminate it from the read name */
        if( NULL != fwd_slash_pntr )
            *fwd_slash_pntr = '\0';
    }
    
    /** If necessary, Determine the read end **/
    if( end == UNKNOWN )
    {
        /* find the final / */
        char* fwd_slash_pntr = 
            strrchr ( readname, '/' );

        /* if there is no slash, assume its not a paired end read */
        if( fwd_slash_pntr == NULL )
        {
            end = NORMAL;
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
            
        }
    }
    
    /* calculate the length of the readname */
    int readname_len = strlen( readname );
    
    /*** get the actual read */
    return_code = fscanf( input_file, "%s\n", read );
    int read_len = strlen( read );

    /*** get the next read name - we discard this */
    /* We can't use  
       return_code = fscanf( input_file, "+%*s\n" );       
       because, if the line is just +\n, it wont find a string and will
       skip to the next line. Therefore, we use this messy while loop.
       TODO - make the assert a real error.
    */

    assert( '+' == getc( input_file)  );
    while( '\n' != getc( input_file  ) );
    
    /*** get the quality score */
    return_code = fscanf( input_file, "%s\n", quality );

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
filter_rawread( struct rawread* r )
{
    /***************************************************************
     * check to make sure this read is mappable
     * we consider a read 'mappable' if:
     * 1) The penalties array is not too low
     *
     */

    /* Make sure the global option has been set 
       ( it's initialized to -1 ); */
    assert( min_num_hq_bps >= 0 );

    /* 
     * count the number of a's and t's. 
     */
    int num_hq_bps = 0;
    int i;
    for( i = 0; i < r->length; i++ )
    {
        double qual = 1 - mutation_probability( 
            ( (unsigned char) (r->error_str)[i] )
        );
        
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

    (*rdb)->lock = malloc( sizeof(pthread_mutex_t) );
    pthread_mutex_init( (*rdb)->lock, NULL );

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
    
    pthread_mutex_destroy( rdb->lock );
    free( rdb->lock );
    
    free( rdb );
    return;
}


void
add_single_end_reads_to_rawread_db( 
    struct rawread_db_t* rdb, char* rifname, 
    enum inputfile_type iftype )
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
    
    return;
}

void
add_paired_end_reads_to_rawread_db( 
    struct rawread_db_t* rdb,
    char* rifname1, char* rifname2,
    enum inputfile_type iftype )
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
get_next_mappable_read_from_rawread_db( 
    struct rawread_db_t* rdb, long* readkey,
    struct rawread** r1, struct rawread** r2 )
{
    int rv = -10;
    
    rv = get_next_read_from_rawread_db( 
            rdb, readkey, r1, r2 );
    
    /* While this rawread is unmappable, continue */ 
    /* 
       This logic is complicated because it has to deal with both paired end
       and non paired end reads. But, basically, it says that r1 cant be null
       and r2 is null and r1 is mappable or r2 is not null ( ie this is paired
       end ) and both r1 and r2 are mappable.
    */
    while( *r1 != NULL 
           && (      ( *r2 == NULL && filter_rawread( *r1 ) == true )
                  || ( filter_rawread( *r1 ) == true 
                       && filter_rawread( *r2 ) == true )
              )
           )
    {
        free_rawread( *r1 );
        if ( *r2 != NULL )
            free_rawread( *r2 );
        
        rv = get_next_read_from_rawread_db( 
                rdb, readkey, r1, r2 );
    }
    
    return rv;
}

int
get_next_read_from_rawread_db( 
    struct rawread_db_t* rdb, long* readkey, 
    struct rawread** r1, struct rawread** r2 )
{
    pthread_mutex_lock( rdb->lock );

    /* 
     * Store the return value - 
     * 0 indicates success, negative failure 
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
            pthread_mutex_unlock( rdb->lock );
            return EOF;
        }        
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
            /* This happends inside the function */
            assert( rawread_db_is_empty( rdb ) );
            *r1 = NULL;
            *r2 = NULL;
            pthread_mutex_unlock( rdb->lock );
            return EOF;
        }
        
        rv = populate_read_from_fastq_file( 
            rdb->paired_end_2_reads, r2, SECOND );
        
        /* if returned an EOF, then it should have been returned for r1 */
        assert( rv == 0 );
    }

    /* increment the read counter */
    *readkey = rdb->readkey;
    rdb->readkey += 1;

    pthread_mutex_unlock( rdb->lock );
    return 0;
}

/******* END raw read DB code ********************************************/
