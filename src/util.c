#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h> /* gettimeofday() */

#include <string.h>
#include <limits.h>

#include <pthread.h>

/* to find out the  number of available processors */
#include <sys/sysinfo.h>
/* mkdir */
#include <sys/stat.h>
/* chdir */
#include <sys/unistd.h>
/* check for system signals */
#include <signal.h>
#include <sys/wait.h>

#include "config.h"
#include "util.h"

#include "log.h"

void
safe_mkdir(char* dir)
{
    int error = mkdir( dir, S_IRWXU | S_IRWXG | S_IRWXO );
    if( -1 == error ) {
        statmap_log( LOG_FATAL, "Cannot mkdir %s", dir );
    }
}

void
safe_copy_into_output_directory( char* fname, char* output_dir, char* output_fname )
{
    int error;

    char* realpath_buffer = realpath( output_dir, NULL );
    assert( NULL != realpath_buffer );

    char buffer[500];
    sprintf( buffer, "cp %s %s/%s", fname, realpath_buffer, output_fname );
    statmap_log( LOG_NOTICE, "Copying '%s' to the output directory", fname );
    error = system( buffer );
    if (WIFSIGNALED(error) &&
        (WTERMSIG(error) == SIGINT || WTERMSIG(error) == SIGQUIT))
    {
        statmap_log( LOG_FATAL, "Copy failure; system call '%s' failed", buffer );
    }
    
    free( realpath_buffer );
    
    return;
}

void
safe_link_into_output_directory( char* fname, char* output_dir, char* output_fname )
{
    int error;

    char* realpath_buffer = realpath( output_dir, NULL );
    assert( NULL != realpath_buffer );

    char* input_realpath_buffer = realpath( fname, NULL );
    assert( NULL != realpath_buffer );

    char buffer[500];
    sprintf( buffer, "ln -s %s %s/%s", input_realpath_buffer, realpath_buffer, output_fname );
    statmap_log( LOG_NOTICE, "Linking '%s' to the output directory", fname );
    error = system( buffer );
    if (WIFSIGNALED(error) &&
        (WTERMSIG(error) == SIGINT || WTERMSIG(error) == SIGQUIT))
    {
        statmap_log( LOG_FATAL, "Link failure; system call '%s' failed", buffer );
    }
    
    free( realpath_buffer );
    free( input_realpath_buffer );
    
    return;
}

FILE* 
open_check_error( char* fname, char* file_mode )
{
    FILE* tmp;
    tmp = fopen( fname, file_mode );
    if( tmp == NULL )
    {
        statmap_log( LOG_FATAL, "Error opening '%s'", fname );
    }
    return tmp;
}

void*
safe_malloc( size_t num_bytes )
{
    assert( num_bytes > 0 );

    void* ma_ptr = malloc( num_bytes );

    if( ma_ptr == NULL )
    {
        statmap_log( LOG_FATAL, "Error malloc'ing %zi bytes",  num_bytes );
        assert( false );
        exit(-1);
    }

    return ma_ptr;
}

/*
   Wrapper for realloc that checks return value
*/
void*
safe_realloc( void* ptr, size_t size )
{
    assert( size > 0 );

    void* ra_ptr;
    ra_ptr = realloc( ptr, size );

    /* man 3 realloc, line 44 - may return NULL if size == 0 */
    if( ra_ptr == NULL && size != 0 )
    {
        statmap_log( LOG_FATAL, "Error realloc'ing %zi bytes.", size );
        assert( false );
        exit(-1);
    }

    return ra_ptr;
}

enum bool
file_is_empty( FILE* fp )
{
    fseek( fp, 0, SEEK_END );
    if( ftell( fp ) == 0)
        return true;

    rewind( fp );
    return false;
}

int
cmp_ints( const int* i1, const int* i2 )
{
    return *i1 - *i2;
}

int
cmp_floats( const float *f1, const float *f2 )
{
    if( (*f1 - *f2) < 0 )
        return -1;

    if( (*f1 - *f2) > 0 )
        return 1;

    return 0;
}
