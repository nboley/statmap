#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <pthread.h>
#include <string.h>

#include <time.h>
#include <sys/time.h> /* gettimeofday() */

/* check for system signals */
#include <signal.h>
#include <sys/wait.h>

#include <unistd.h>
#include <sys/syscall.h> /* gettid */

#include "config.h"
#include "log.h"
#include "util.h"

/* Statically initialize the mutex when it is declared */
static pthread_mutex_t log_mutex = PTHREAD_MUTEX_INITIALIZER;
static FILE* log_fp = NULL;
static char* tmp_logfile_path = NULL;

/* all log messages with level >= this level should also be printed to stderr
 * so the user can respond appropriately */
enum LOG_LEVEL nontrivial_log_level;

void
init_initial_logging()
{
    /* We need somewhere to write logging output until the "real" log file is
     * created with init_logging (which cannot happen until the Statmap output
     * directory has been created. */
    
    /* Build a unique filename with the current time and a random nonce (to
     * avoid the case where multiple Statmap processes are started at the same
     * time, say, by a script, which could cause them all to overwrite each
     * other's temporary log files) */
    time_t rawtime = time( NULL );
    struct tm* timeinfo = localtime( &rawtime );

    char buffer[255];
    char* buf_pos = buffer;
    buf_pos += strftime( buffer, 255, "/tmp/statmap-%Y_%m_%d_%H_%M_%S",
            timeinfo );
    buf_pos += sprintf( buf_pos, "-%i.tmp.log", rand() );
    tmp_logfile_path = calloc( strlen(buffer)+1, sizeof(char) );
    strncpy( tmp_logfile_path, buffer, strlen(buffer) );

    assert( log_fp == NULL );
    log_fp = fopen( tmp_logfile_path, "w" );
    if( log_fp == NULL )
    {
        /* log directly to stderr, since we couldn't open the logfile */
        fprintf( stderr, "FATAL : Could not open temporary log file at '%s'\n",
                tmp_logfile_path );
        assert(false);
        exit(1);
    }        
}

void
init_logging( enum LOG_LEVEL min_nontrivial_log_level )
{
    int rv;

    /* Check to see if the log_fp is already initialized to the temporary fp.
     * If it is, close it, then copy its contents to the new logfile and delete
     * the temporary file. If this code is being called by a utility, it might
     * not have been initialized. */
    if( log_fp != NULL )
    {
        rv = fclose( log_fp );
        assert( rv == 0 );

        /* Do the copy here so we know how to handle logging (exceptional,
         * since no log file is open at the moment) */
        char cmd[500];
        sprintf( cmd, "cp %s %s", tmp_logfile_path, LOG_FNAME );
        int error;
        error = system( cmd );
        if (WIFSIGNALED(error) &&
                (WTERMSIG(error) == SIGINT || WTERMSIG(error) == SIGQUIT))
        {
            fprintf( stderr, 
                     "FATAL : Could not copy temporary log file '%s' to log file '%s'\n",
                     tmp_logfile_path, LOG_FNAME );
            assert(false);
            exit(1);
        }

        free( tmp_logfile_path );
    }

    /* Initialize the log fp. We append to the file because it already exists
     * (either copied from the temporary log file, or possibly created by an
     * initial Statmap run and is now being logged to by a utility using
     * functions in the shared library */

    /* Assume we are already ch'dired into a Statmap output directory */
    log_fp = fopen( LOG_FNAME, "a" );
    if( log_fp == NULL )
    {
        /* log directly to stderr, since we couldn't open the log file */
        fprintf( stderr, "FATAL : Could not open log file at '%s'\n",
                 LOG_FNAME );
        assert(false);
        exit(1);
    }
    
    nontrivial_log_level = min_nontrivial_log_level;

}

static inline const char*
log_level_to_string( enum LOG_LEVEL log_level )
{
    switch(log_level)
    {
        case LOG_DEBUG: return "DEBUG";
        case LOG_INFO:  return "INFO";
        case LOG_NOTICE: return "NOTICE";
        case LOG_WARNING: return "WARNING";
        case LOG_ERROR: return "ERROR";
        case LOG_FATAL: return "FATAL";
        default:
            assert(0);
            return "UNKNOWN";
    }
}

void
statmap_log( enum LOG_LEVEL log_level, const char* format, ... )
{
    /* Build timestamp string */
    char timestamp[100];
    time_t rawtime = time( NULL );
    struct tm* timeinfo;
    if( (time_t)-1 != rawtime ) {
        timeinfo = localtime( &rawtime );
        strftime(timestamp, 100, "%Y-%m-%d %H:%M:%S", timeinfo );
    } else {
        strcpy( timestamp, "UNKNOWN-TIME" );
    }

    /* Build the log message format string */
    char log_msg[1024];
    va_list args;
    va_start(args, format);
    vsprintf(log_msg, format, args);
    va_end(args);

    /* Get the thread id from the kernel */
    pid_t thread_id = (pid_t) syscall(SYS_gettid);

    /* Build the log line to output (prepends thread_id, log level) */
    char log_line[1024];
    sprintf( log_line, "[%s] [%d] %-10s: %s\n", timestamp, thread_id,
             log_level_to_string(log_level), log_msg );

    /* Print the message to the log file, if one is open */
    if( NULL != log_fp  )
    {
        /* acquire a mutex to make sure that the log isnt garbled */
        pthread_mutex_lock( &log_mutex );
        fputs( log_line, log_fp );
        pthread_mutex_unlock( &log_mutex );
        fflush( log_fp );
    }
    
    /* Print to stderr if there is no log file or the message is non-trivial */
    if( log_fp == NULL || log_level >= nontrivial_log_level )
    {
        fputs( log_line, stderr );
    }    
    
    /* If this was a fatal log message, exit */
    if( log_level == LOG_FATAL )
    {
        assert( false );
        exit(1);
    }
}

void
finish_logging()
{
    fclose( log_fp );
}
