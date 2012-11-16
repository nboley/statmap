#include <stdio.h>
#include <stdarg.h>
#include <pthread.h>

#ifndef LOG_H
#define LOG_H

enum LOG_LEVEL {
    LOG_DEBUG = 1,
    LOG_INFO = 2,
    LOG_NOTICE = 3,
    LOG_WARNING = 4,
    LOG_ERROR = 5,
    LOG_FATAL = 6
};

#define LOG_FNAME "log"
/* all log messages with level >= this level should also be printed to stderr
 * so the user can respond appropriately */
#define NONTRIVIAL_LOG_LEVEL LOG_INFO

void
init_initial_logging();

void
init_logging();

void
statmap_log( enum LOG_LEVEL log_level, const char* format, ... );

void
finish_logging();

#endif // LOG_H
