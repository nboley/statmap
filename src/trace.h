#ifndef TRACE_H
#define TRACE_H 

#define TRACE_TYPE float

#define TRACE_MAGIC_NUMBER "BTRACE"

#include "genome.h"

// define this to use mutexes, otherwise use spinlocks
// #define USE_MUTEX
#define USE_SPINLOCK

struct trace_segment_t {
    int real_track_id;
    int real_chr_id;
    int real_start;

    int length;
    TRACE_TYPE* data;
    pthread_mutex_t* data_lock;
};

struct trace_segments_t {
    int num_segments;
    struct trace_segment_t* segments;
};

struct trace_t {
    int num_tracks;
    char** track_names;
    
    int num_chrs;
    char** chr_names;
    unsigned int* chr_lengths;

    // num_tracks x num_chrs   
    struct trace_segments_t** segments;
};

void
init_trace( struct genome_data* genome,
            struct trace_t** traces,
            const int num_traces,
            char** track_names );

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
write_wiggle_from_trace_to_stdout(
    struct trace_t* traces,
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

/********************************************************************
 * Trace aggregate functions
 *
 ********************************************************************/

TRACE_TYPE
trace_agg_sum( const TRACE_TYPE a, const TRACE_TYPE b );

TRACE_TYPE
trace_agg_min( const TRACE_TYPE a, const TRACE_TYPE b );

TRACE_TYPE
trace_agg_max( const TRACE_TYPE a, const TRACE_TYPE b );

void
update_trace_segments_from_mapped_read_array(
    struct trace_segments_t* trace_segments,
    float* update_vals,
    float scale_factor,
    int start,
    int stop
);

/* Wrapper function to update from a uniform kernel */
void
update_trace_segments_from_uniform_kernel(
    struct trace_segments_t* trace_segments,
    float scale_factor,
    int start,
    int stop
);

double
accumulate_from_trace(
    const struct trace_t* const traces,
    int track_index,
    int chr_index,
    int start,
    int stop
);

double
accumulate_from_traces(
    const struct trace_t* const traces,
    int chr_index,
    int start,
    int stop
);

#endif // #define TRACE_H
