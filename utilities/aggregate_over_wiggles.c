#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "../src/wiggle.h"

int min_num_hq_bps = -1;
int num_threads = 1;

void usage()
{
    fprintf(stderr, "Usage: ./aggregate_over_wiggles (min|max|sum) file(s).wig > output.wig\n");
}


int 
main( int argc, char** argv )
{
    if( argc <= 2 )
    {
        usage();
        exit( -1 );
    }

    /* determine the aggregate type */
    float (*agg_fn)( const struct wig_line_info*, const int, const int ) = NULL;

    if( 0 == strcmp("max", argv[1] ) )
    {
        agg_fn = wig_lines_max;
    } else if( 0 == strcmp("min", argv[1] ) ) {
        agg_fn = wig_lines_min;
    } else if( 0 == strcmp("sum", argv[1] ) ) {
        agg_fn = wig_lines_sum;
    } else {
        fprintf( stderr, "FATAL     : Unrecognized aggregate '%s'\n", argv[1] );
        exit( -1 );
    }

    /* open the wiggles */
    FILE** wigs = calloc(argc-2, sizeof(FILE*));
    int i;
    for( i = 0; i < argc-2; i++ )
    {
        wigs[i] = fopen( argv[i+2], "r" );
        if( NULL == wigs[i] )
        {
            perror( "FATAL     : Could not open wiggle file for reading" );
            exit( -1 );
        }
    }

    aggregate_over_wiggles( wigs, argc-2, stdout, FLT_EPSILON, agg_fn );
    
    for( i = 0; i < argc-2; i++ )
        fclose( wigs[i] );
    
    free( wigs );

    return 0;
}


