#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <pthread.h>

#include "config.h"
#include "log.h"

/* Statically initialize the mutex when it is declared */
static pthread_mutex_t log_mutex = PTHREAD_MUTEX_INITIALIZER;
static FILE* log_fp;

void
init_logging()
{
    /* Initialize the log fp. We append to the file in case it already exists
     * (for example, if it was generated during marginal mapping and we run
     * a utility that does iterative mapping) */

    /* ATM, we assume we have already ch'dired into the Statmap output
     * directory */
    //log_fp = fopen( LOG_FNAME, "a" );
    log_fp = stderr;
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

/* TODO - add fn name, thread id to log output
 * gettid/getpid for thread id. If program is single threaded, tid == pid. */

void
statmap_log( enum LOG_LEVEL log_level, char* fmt, ... )
{
    /* Build the log message format string */
    char log_msg[1024];
    va_list args;
    va_start(args, fmt);
    vsprintf(log_msg, fmt, args);
    va_end(args);

    /* Build the log line to output (prepends log level) */
    char log_line[1024];
    sprintf( log_line, "%-10s: %s\n",
             log_level_to_string(log_level), log_msg );

    /* Print the message to the log file (and stderr if non-trivial) */
    pthread_mutex_lock( &log_mutex );
    fputs( log_line, log_fp );
    /*
     * TODO - come back to this once we've properly set up logging to a file
    if( log_level >= NONTRIVIAL_LOG_LEVEL )
    {
        fputs( log_line, stderr );
    }
    */
    pthread_mutex_unlock( &log_mutex );

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
