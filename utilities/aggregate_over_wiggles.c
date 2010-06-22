
#include <stdio.h>
#include <stdlib.h>

#include "../src/wiggle.h"

void usage()
{
    fprintf(stderr, "Usage: ./aggregate_over_wiggles files.wig > output.wig\n");
}


int 
main( int argc, char** argv )
{
    if( argc <= 1 )
    {
        usage();
        exit( -1 );
    }

    /* open the wiggles */
    FILE** wigs = calloc(argc-1, sizeof(FILE*));
    int i;
    for( i = 1; i < argc; i++ )
    {
        wigs[i] = fopen( argv[i], "r" );
        if( NULL == wigs[i] )
        {
            perror( "FATAL     : Could not open wiggle file for reading" );
            exit( -1 );
        }
    }

    aggregate_over_wiggles( wigs, argc-1, stdout );
    
    for( i = 1; i < argc; i++ )
        fclose( wigs[i] );
    
    free( wigs );
}


