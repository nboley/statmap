#ifndef TRACE_H
#define TRACE_H

#define TRACE_TYPE float

#define TRACE_MAGIC_NUMBER "BTRACE"

/* Trace spinlock granularity */
/* How many basepairs are grouped together for a single lock */
#define TM_GRAN 5000

struct genome_data;

// define this to use mutexes, otherwise use spinlocks
// #define USE_MUTEX
#define USE_SPINLOCK

struct trace_t {
    int num_tracks;
    char** track_names;
    
    int num_chrs;
    char** chr_names;
    unsigned int* chr_lengths;
    
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
init_trace_array( size_t size );

void
init_trace( struct genome_data* genome,
            struct trace_t** traces,
            const int num_traces,
            char** track_names );

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
multiply_trace_by_scalar( struct trace_t* traces, double value );

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
aggregate_over_trace_pairs(  
    struct trace_t* update_trace, 
    const struct trace_t* const other_trace,
    TRACE_TYPE (*aggregate)( const TRACE_TYPE, const TRACE_TYPE )
);

void
aggregate_over_traces_from_fnames(  
    char** fnames,
    int num_files,
    struct trace_t** agg_trace,
    TRACE_TYPE (*aggregate)( const TRACE_TYPE, const TRACE_TYPE )
);

/****************************************************************
 * Wiggle IO
 *
 ****************************************************************/

extern void
write_wiggle_from_trace_to_stream( 
    struct trace_t* traces,
    FILE* os, /* output stream */                           
    const double filter_threshold 
);

extern void
write_wiggle_from_trace( 
    struct trace_t* traces,
    const char* output_fname,                           
    const double filter_threshold 
);


/****************************************************************
 * Pickling and Unpickling 
 *
 ****************************************************************/

void
write_trace_to_stream( struct trace_t* trace, FILE* os );

void
write_trace_to_file( struct trace_t* trace, char* fname );

void
load_trace_from_stream( struct trace_t** trace, FILE* is );

void
load_trace_from_file( struct trace_t** trace, char* fname );

#endif // #define TRACE_H
