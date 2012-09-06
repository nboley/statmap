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

void
safe_mkdir(char* dir)
{
    int error = mkdir( dir, S_IRWXU | S_IRWXG | S_IRWXO );
    if( -1 == error )
    {
        fprintf(stderr, "FATAL       :  Cannot make %s\n", dir );
        exit( -1 );
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
    fprintf(stderr, "NOTICE      :  Copying '%s' to the output directory\n",  fname );
    error = system( buffer );
    if (WIFSIGNALED(error) &&
        (WTERMSIG(error) == SIGINT || WTERMSIG(error) == SIGQUIT))
    {
        fprintf(stderr, "FATAL     : Failed to call '%s'\n", buffer );
        perror( "Copy failure");
        exit( -1 );
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
    fprintf(stderr, "NOTICE      :  Linking '%s' to the output directory\n",  fname );
    error = system( buffer );
    if (WIFSIGNALED(error) &&
        (WTERMSIG(error) == SIGINT || WTERMSIG(error) == SIGQUIT))
    {
        fprintf(stderr, "FATAL     : Failed to call '%s'\n", buffer );
        perror( "Link failure");
        exit( -1 );
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
        fprintf( stderr, "Error opening '%s\n'", fname );
        exit( -1 );
    }
    return tmp;
}

/*
   Wrapper for realloc that checks return value
*/
void*
safe_realloc( void* ptr, size_t size )
{
    void* ra_ptr;
    ra_ptr = realloc( ptr, size );

    /* man 3 realloc, line 44 - may return NULL if size == 0 */
    if( ra_ptr == NULL && size != 0 )
    {
        fprintf(stderr, "FATAL       :  Error realloc'ing %zi bytes.\n", size);
        exit( -1 );
    }

    return ra_ptr;
}

enum bool
file_is_empty( FILE* fp )
{
    fseek( fp, 0, SEEK_END );
    if( ftell( fp ) == 0)
        return true;

    return false;
}
