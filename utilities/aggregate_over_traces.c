#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <assert.h>

#include "../src/trace.h"

void usage()
{
    fprintf(stderr, "Usage: ./aggregate_over_traces (min|max|sum) output.bin.trace file(s).bin.trace\n");
}

int 
main( int argc, char** argv )
{
    if( argc <= 3 )
    {
        usage();
        exit( -1 );
    }

    /* determine the aggregate type */
    TRACE_TYPE (*agg_fn)( const TRACE_TYPE, const TRACE_TYPE ) = NULL;

    if( 0 == strcmp("max", argv[1] ) )
    {
        agg_fn = trace_max_agg;
    } else if( 0 == strcmp("min", argv[1] ) ) {
        agg_fn = trace_min_agg;
    } else if( 0 == strcmp("sum", argv[1] ) ) {
        agg_fn = trace_sum_agg;
    } else {
        fprintf( stderr, "FATAL     : Unrecognized aggregate '%s'\n", argv[1] );
        exit( -1 );
    }

    /* load the first trace as the trace that we will update */
    struct trace_t* update_trace;
    load_trace_from_file( &update_trace, argv[3] );
    
    int i;
    for( i = 1; i < argc-3; i++ )
    {
        struct trace_t* curr_trace;
        load_trace_from_file( &curr_trace, argv[i+3] );
        aggregate_over_trace_pairs( update_trace, curr_trace, agg_fn );
        close_traces( curr_trace );
    }
    
    write_trace_to_file( update_trace, argv[2] );
    
    close_traces( update_trace );
    
    return 0;
}

