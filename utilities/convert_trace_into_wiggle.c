#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <assert.h>

#include "../src/trace.h"

void usage()
{
    fprintf(stderr, "Usage: ./convert_trace_into_wiggle input.bin.trace [filter_threshold] > output.wig\n");
}

int 
main( int argc, char** argv )
{
    if( argc <= 1 || argc >= 4 )
    {
        usage();
        exit( -1 );
    }

    /* load the trace into memory */
    struct trace_t* trace;
    load_trace_from_file( &trace, argv[1] );

    double filter_threshold = 0;
    if( argc == 3 )
    {
        filter_threshold = atof( argv[2] );
    }
    
    write_wiggle_from_trace_to_stream( trace, stdout, filter_threshold );
    
    close_traces( trace );
    
    return 0;
}


