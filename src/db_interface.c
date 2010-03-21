/* Copyright (c) 2009-2010, Nathan Boley */

#include <sys/types.h>
#include <sys/stat.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>

#include "config.h"
#include "statmap.h"
#include "db_interface.h"


/*** Generic Candidate Mappings DB Code *********************************************/

inline FILE* open_check_error( char* fname, char* file_mode )
{
    FILE* tmp;
    tmp = fopen( fname, file_mode );
    if( tmp == NULL )
    {
        fprintf( stderr, "Error opening '%s'", fname );
        exit( -1 );
    }
    return tmp;
}

/* TODO - think about this */

int
cmp_joining_queue_datum( const void* a, const void* b )
{
    /* a is of type joining_queue_datum** */
    joining_queue_datum** c = (joining_queue_datum**) a;
    joining_queue_datum** d = (joining_queue_datum**) b;

    /* Make NULL's compare largest */
    if( *c == NULL && *d == NULL )
        return 0;
    if( *d == NULL )
        return -1;
    if( *c == NULL )
        return 1;
    
    if( (*c)->read_id != (*d)->read_id )
        return (*c)->read_id - (*d)->read_id;
    
    return cmp_candidate_mappings( &((*c)->mapping), 
                                   &((*d)->mapping)  );

}

void
init_candidate_mappings_db( candidate_mappings_db* db,
                            char* candidate_mappings_prefix)
{
    int error;
    
    char* cmdir_name;
    if( candidate_mappings_prefix == NULL )
    {
        char dir_name_template[100] = "cand_mappings_XXXXXX";
        cmdir_name =  mkdtemp(dir_name_template);
        if( cmdir_name == NULL )
        {
            fprintf( stderr, 
                     "FATAL       :  Could not create temp directory with template '%s'\n", 
                     dir_name_template );
            exit( -1 );
        }            
    } else {
        cmdir_name = candidate_mappings_prefix;
        int rv = mkdir( cmdir_name, 700 );
        if( rv != 0 )
        {
            fprintf( stderr,
                     "FATAL       :  Failed (errno %i) to create dir '%s'\n", 
                     rv, cmdir_name );
            exit( -1 );
        }
    }
    
    fprintf( stderr, 
             "NOTICE      :  Created Directory for Candidate Mappings: %s\n", 
             cmdir_name 
    );
    db->data_dir = malloc( sizeof(char)*(strlen( cmdir_name )+ 1) );
    memcpy( db->data_dir, cmdir_name, sizeof(char)*(strlen( cmdir_name )+ 1) );

    char tfname[1000];

    /* If we want, initialize BDB */
    db->bdb = NULL;
    
    /* TODO - allow for postgres */
    // db->pdb = NULL;
    
    /* TODO - allow for buffered reads */
    db->buffer = NULL;

    db->full_fwd_rds = malloc(num_threads*sizeof(FILE*));
    db->full_bkwd_rds = malloc(num_threads*sizeof(FILE*));
    db->full_1st_fwd_rds = malloc(num_threads*sizeof(FILE*));
    db->full_1st_bkwd_rds = malloc(num_threads*sizeof(FILE*));
    db->full_2nd_fwd_rds = malloc(num_threads*sizeof(FILE*));
    db->full_2nd_bkwd_rds = malloc(num_threads*sizeof(FILE*));
    
    error = chdir( cmdir_name );
    assert( error == 0 );
    
    int i = 0;
    for( i = 0; i < num_threads; i++ )
    {
        sprintf( tfname, "thread%i.full_fwd_reads.mapped",
                 i);
        db->full_fwd_rds[i] = open_check_error( tfname, "w+" );

        sprintf( tfname, "thread%i.full_bkwd_reads.mapped",
                 i);
        db->full_bkwd_rds[i] = open_check_error( tfname, "w+" );

        sprintf( tfname, "thread%i.full_1st_fwd_reads.mapped",
                 i);
        db->full_1st_fwd_rds[i] = open_check_error( tfname, "w+" );

        sprintf( tfname, "thread%i.full_1st_bkwd_reads.mapped",
                 i);
        db->full_1st_bkwd_rds[i] = open_check_error( tfname, "w+" );

        sprintf( tfname, "thread%i.full_2nd_fwd_reads.mapped",
                 i);
        db->full_2nd_fwd_rds[i] = open_check_error( tfname, "w+" );

        sprintf( tfname, "thread%i.full_2nd_bkwd_reads.mapped",
                 i);
        db->full_2nd_bkwd_rds[i] = open_check_error( tfname, "w+" );
    }

    /* TODO - get curr dir and change back to that */
    error = chdir( "../" );
    assert( error == 0 );
 
    return;
}

void
close_candidate_mappings_db( candidate_mappings_db* db )
{  
    int error;

    int i;
    for( i = 0; i < num_threads; i++ )
    {
        fclose( db->full_fwd_rds[i] );
        fclose( db->full_bkwd_rds[i] );
        fclose( db->full_1st_fwd_rds[i] );
        fclose( db->full_1st_bkwd_rds[i] );
        fclose( db->full_2nd_fwd_rds[i] );
        fclose( db->full_2nd_bkwd_rds[i] );
    }

    free( db->full_fwd_rds );
    free( db->full_bkwd_rds );
    free( db->full_1st_fwd_rds );
    free( db->full_1st_bkwd_rds );
    free( db->full_2nd_fwd_rds );
    free( db->full_2nd_bkwd_rds );

    /* Remove the cand mapping files */
    char buffer[255];
    sprintf( buffer, "rm %s -rf\n", db->data_dir );
    fprintf( stderr, "NOTICE      :  Removing tmp directoty: %s", buffer );
    error = system( buffer );
    if( error != 0 )
    {
        fprintf(stderr, "FATAL       :  Error removing temp files.\n");
        fprintf(stderr, "Failed cmd: %s\n", buffer);
        assert(false);
        exit( -1 );
    }


    free( db->data_dir );
    
    return;
}


inline FILE*
determine_fp_from_candidate_mapping( 
    candidate_mapping* mapping, 
    candidate_mappings_db* db,
    int thread_num )
{
    if( mapping->rd_type == SINGLE_END )
    {
        if( mapping->rd_strnd == FWD ) {
            return db->full_fwd_rds[thread_num];
        } else {
            assert( mapping->rd_strnd == BKWD );
            return db->full_bkwd_rds[thread_num];
        }
    } else 
    if( mapping->rd_type == PAIRED_END_1 )
    {
        if( mapping->rd_strnd == FWD ) {
            return db->full_1st_fwd_rds[thread_num];
        } else {
            assert( mapping->rd_strnd == BKWD );
            return db->full_1st_bkwd_rds[thread_num];
        }
    } else 
    {   assert( mapping->rd_type == PAIRED_END_2 );
            
        if( mapping->rd_strnd == FWD ) {
            return db->full_2nd_fwd_rds[thread_num];
        } else {
            assert( mapping->rd_strnd == BKWD );
            return db->full_2nd_bkwd_rds[thread_num];
        }
    }
}

void
add_candidate_mappings_to_db ( 
    candidate_mappings_db* db, 
    candidate_mappings* mappings,
    int thread_num,
    char* readname,
    unsigned long read_id )
{
    if ( mappings->length == 0 )
        return;

    int error;
    
    /* TODO - add proper DB support */
    /* Determine the correct file to write this read to  */        
    /* BUG - FIXME */
    FILE* fp =  
        determine_fp_from_candidate_mapping( 
            mappings->mappings, db, thread_num );
        
    /* write the key to the file */
    error = fprintf( fp, "%s\t%lu", readname, read_id );
    error = fputc( '\0', fp );
    
    unsigned int i;
    for( i = 0; i < mappings->length; i++ )
    {
        /* Move to the correct candidate mapping */
        candidate_mapping* mapping = mappings->mappings + i;

        error = putc( '\0', fp );
        error = fwrite( mapping, sizeof(candidate_mapping), 1, fp );
    }
    error = putc( '\n', fp );
    
    return;
}

void
close_candidate_mappings_cursor( 
    candidate_mappings_db_cursor* cursor )
{
    int i;
    for( i = 0; i < cursor->queue_len; i++ )
    {
        /* free the queue_datum ) */
        if( (cursor->mappings_queue)[i] != NULL )
        {
            free( (cursor->mappings_queue)[i] );
        }
    }
    
    free( cursor->mappings_queue );

    /* free the cursor memory */
    free( cursor );
}

inline int
get_readname_from_stream( char* buffer, 
                          unsigned long* read_id, 
                          int max_chars, 
                          FILE* stream )
{
    int error;
    int i;
    for( i = 0; i < max_chars-1; i++ )
    {
        buffer[i] = fgetc( stream );
        /* If we are at the end of the stream */
        if( feof( stream ) )
        {
            buffer[i] = '\0';
            return EOF;
        }

        /* If we are at the end of the readname */
        if( buffer[i] == '\t' )
        {
            /* Set the NULL byte */
            buffer[i] = '\0';
            
            /* Read in the read id */
            error = fscanf( stream, "%lu", read_id );
            assert( error == 1 );
            fgetc( stream );
            
            return 0;
        }
    }
    buffer[ max_chars - 1 ] = '\0';

    fprintf(stderr, "Failed to read readname in 'get_readname_from_stream'\n" );
    exit( -1 );
}

void
open_candidate_mappings_cursor (
    candidate_mappings_db* db,
    candidate_mappings_db_cursor** cursor)
{
    int i;
    int error;
    
    *cursor = malloc( sizeof( candidate_mappings_db_cursor ) );
    (*cursor)->db = db;
    (*cursor)->queue_len = 6*num_threads;
    (*cursor)->mappings_queue 
        = malloc(sizeof(joining_queue_datum*)*(*cursor)->queue_len);
    
    FILE** fps = malloc(sizeof(FILE*)*6*num_threads);
    for( i = 0; i < num_threads; i++ )
    {
        fps[6*i + 0] = db->full_fwd_rds[i];
        fps[6*i + 1] = db->full_bkwd_rds[i];
        fps[6*i + 2] = db->full_1st_fwd_rds[i];
        fps[6*i + 3] = db->full_1st_bkwd_rds[i];
        fps[6*i + 4] = db->full_2nd_fwd_rds[i];
        fps[6*i + 5] = db->full_2nd_bkwd_rds[i];
    };

    /* initialize the queue */    
    for( i = 0; i < (*cursor)->queue_len; i++ )
    {
        /* Move the fp's to the beginning of the stream */ 
        rewind( fps[i] );

        /* initialize space for the queue */
        ((*cursor)->mappings_queue)[i] 
            = malloc( sizeof( joining_queue_datum ) );

        /* add the file to the stream */
        ((*cursor)->mappings_queue)[i]->stream = fps[i];
        
        /* Get the readnames from the stream while the reads 
           dont have mappings */
        do { 
            error = get_readname_from_stream( 
                ((*cursor)->mappings_queue)[i]->key,
                &(((*cursor)->mappings_queue)[i]->read_id),
                MAX_KEY_SIZE + 1,
                fps[i]
            );
            
            error = getc( fps[i] );
        } while( error == '\n' );
        
        /* if the file is at eof, set the queue to NULL */
        if( EOF == error )
        {
            free( ((*cursor)->mappings_queue)[i] );
            ((*cursor)->mappings_queue)[i] = NULL;
          
        } else {        
            /* read in the mapped location and add it into the queue */
            error = fread( &((((*cursor)->mappings_queue)[i])->mapping), 
                    sizeof(candidate_mapping), 1, fps[i] );
            
        }
    }

    /* sort the queue */
    qsort( (*cursor)->mappings_queue, 
           (*cursor)->queue_len, 
           sizeof( joining_queue_datum* ), 
           cmp_joining_queue_datum );

    free( fps );
    
    return;
}

int
get_next_candidate_mapping_from_cursor( 
    candidate_mappings_db_cursor* cursor, 
    candidate_mappings** mappings,
    char* key )
{
    /* initialize a place to store the reads */
    init_candidate_mappings( mappings );
    
    /* sort the queue */
    qsort( cursor->mappings_queue, 
           cursor->queue_len, 
           sizeof( joining_queue_datum* ), 
           cmp_joining_queue_datum );
    
    /* check to see if we are done */
    if( (cursor->mappings_queue)[0] == NULL )
    {
        free_candidate_mappings( *mappings );
        return CURSOR_EMPTY;
    }

    /* populate the key with the first key in the queue */
    int key_len = strlen( ((cursor->mappings_queue)[0])->key ) + 1;
    memcpy( key, ((cursor->mappings_queue)[0])->key, key_len );
    /* explicitly set the null byte */
    key[key_len-1] = '\0';
    
    /* now, keep adding data until the key changes */
    do
    {
        /* Add the current mapped location into mapped locations */
        add_candidate_mapping( *mappings, 
                               &(((cursor->mappings_queue)[0])->mapping) );
        
        /* get the next char in the top file to check if we are on a new read */
        /* If we are, this should be a new line */
        char rv = getc( ((cursor->mappings_queue)[0])->stream );
        
        /* If we have moved to the next read */
        if( rv == '\n' )
        {
            /* loop until we have a read with locations, or an EOF */
            /* This is to ensure that an active read is always in the queue */
            while( rv == '\n' )
            {
                rv = get_readname_from_stream( 
                    (cursor->mappings_queue)[0]->key,
                    &((cursor->mappings_queue)[0]->read_id),
                    MAX_KEY_SIZE + 1,
                    (cursor->mappings_queue)[0]->stream
                );
                
                rv = getc( ((cursor->mappings_queue)[0])->stream );
            }
        }
                
        if( rv == EOF )
        {
            assert( feof( ((cursor->mappings_queue)[0])->stream ) );
            free( (cursor->mappings_queue)[0] );
            (cursor->mappings_queue)[0] = NULL;
            // return 0;
        } else {
            /* read in the next mapping */
            assert( rv == '\0' );
            
            rv = fread( &((cursor->mappings_queue)[0]->mapping),
                        sizeof(candidate_mapping), 1,
                        ((cursor->mappings_queue)[0])->stream
            );
        }      
        
        /* finally, resort the queue with the new item */
        qsort( cursor->mappings_queue, 
               cursor->queue_len, 
               sizeof( joining_queue_datum* ), 
               cmp_joining_queue_datum );
        
    } while( (cursor->mappings_queue)[0] != NULL 
             && 0 == strcmp( key, ((cursor->mappings_queue)[0])->key ));

    return 0;
}




