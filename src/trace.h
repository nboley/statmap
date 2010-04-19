#define TRACE_TYPE float

struct genome_data;

typedef struct {
    int num_traces;
    int num_chrs;
    unsigned int* trace_lengths;
    /* num_tracesXnum_chrs matrix */
    TRACE_TYPE*** traces;
} traces_t;

/* build an mmapped array to store the density */
TRACE_TYPE*
init_trace( size_t size );

void
init_traces( struct genome_data* genome,
             traces_t** traces,
             const int num_traces );

void
close_traces( traces_t* traces );

double
sum_traces( traces_t* traces );

void
zero_traces( traces_t* traces );

void
aggregate_over_traces(  traces_t* update_trace, 
                        const traces_t* const other_trace,
                        TRACE_TYPE (*aggregate)( const TRACE_TYPE, const TRACE_TYPE )
    );


void
write_wiggle_from_traces( traces_t* traces,
                          /* the index of the trace to use ( ie, 0 for the fwd trace ) */
                          int trace_index,
                          
                          char** trace_names,
                          
                          const char* output_fname, 
                          const char* track_name,
                          
                          const double filter_threshold );
