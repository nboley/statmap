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

