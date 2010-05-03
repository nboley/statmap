#define TRACE_TYPE float

struct genome_data;

struct trace_t {
    int num_traces;
    int num_chrs;
    unsigned int* trace_lengths;
    /* num_tracesXnum_chrs matrix */
    TRACE_TYPE*** traces;
};

/* build an mmapped array to store the density */
TRACE_TYPE*
init_trace( size_t size );

void
init_traces( struct genome_data* genome,
             struct trace_t** traces,
             const int num_traces );

void
copy_traces( struct trace_t** traces,
             struct trace_t* original );

void
close_traces( struct trace_t* traces );

double
sum_traces( struct trace_t* traces );

void
zero_traces( struct trace_t* traces );

void
set_trace_to_uniform( struct trace_t* traces, double value );

void
aggregate_over_traces(  struct trace_t* update_trace, 
                        const struct trace_t* const other_trace,
                        TRACE_TYPE (*aggregate)( const TRACE_TYPE, const TRACE_TYPE )
    );

extern void
write_wiggle_from_trace( struct trace_t* traces,

                         /* These are usually chr names */
                         char** scaffold_names,
                         /* The names of the various tracks */
                         /* Use null for the indexes */
                         char** track_names,
                         
                         const char* output_fname,                           
                         const double filter_threshold 
    );
