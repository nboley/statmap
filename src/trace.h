#ifndef TRACE_H
#define TRACE_H


#define TRACE_TYPE float

/* Trace spinlock granularity */
/* How many basepairs are grouped together for a single lock */
#define TM_GRAN 5000

struct genome_data;

// define this to use mutexes, otherwise use spinlocks
// #define USE_MUTEX
#define USE_SPINLOCK

struct trace_t {
    int num_traces;
    int num_chrs;
    unsigned int* trace_lengths;
    /* num_traces X num_chrs matrix */
    TRACE_TYPE*** traces;

    #ifdef USE_MUTEX
    pthread_mutex_t*** locks;
    #else
    pthread_spinlock_t*** locks;
    #endif
};

/* build an mmapped array to store the density */
TRACE_TYPE*
init_trace( size_t size );

void
init_traces( struct genome_data* genome,
             struct trace_t** traces,
             const int num_traces );

void
copy_trace_data( struct trace_t* traces,
                 struct trace_t* original );

void
copy_trace( struct trace_t** traces,
            struct trace_t* original );

void
copy_trace_structure( struct trace_t** traces,
                      struct trace_t* original );

void
close_traces( struct trace_t* traces );

void
divide_trace_by_sum( struct trace_t* traces, double value );

void
normalize_traces( struct trace_t* traces );

double
sum_traces( struct trace_t* traces );

void
zero_traces( struct trace_t* traces );

void
set_trace_to_uniform( struct trace_t* traces, double value );

void
apply_to_trace( struct trace_t* traces, double (*fun)(double) );

void
aggregate_over_traces(  struct trace_t* update_trace, 
                        const struct trace_t* const other_trace,
                        TRACE_TYPE (*aggregate)( const TRACE_TYPE, const TRACE_TYPE )
    );

#endif // #define TRACE_H
