#define TRACE_TYPE float

/* Trace spinlock granularity */
/* How many basepairs are grouped together for a single spinlock */
/* Always just call lock_spinlock - this should never be used directly */
#define TM_GRAN 200

struct genome_data;

// define this to use mutexes, otherwise use spinlocks
#define USE_MUTEX
// #define USE_SPINLOCK

struct trace_t {
    int num_traces;
    int num_chrs;
    unsigned int* trace_lengths;
    /* num_tracesXnum_chrs matrix */
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
copy_trace_structure( struct trace_t** traces,
                      struct trace_t* original );

void
close_traces( struct trace_t* traces );

void
normalize_traces( struct trace_t* traces );

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
