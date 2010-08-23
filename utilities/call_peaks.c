#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int num_threads = -1;
int min_num_hq_bps = -1;

#include "../src/wiggle.h"


void usage()
{
    fprintf(stderr, "Usage: ./call_peaks bs_sample_wiggles > called_peaks.wig\n");
}


int 
main( int argc, char** argv )
{
    if( argc < 2 )
    {
        usage();
        exit( -1 );
    }

    /* open the wiggles */
    FILE** ip_wigs = calloc((argc-1)/2, sizeof(FILE*));
    int ip_cntr = -1;
    FILE** nc_wigs = calloc((argc-1)/2, sizeof(FILE*));
    int nc_cntr = -1;
    int i;
    for( i = 0; i < argc-1; i++ )
    {
        int cntr = -1;
        
        /* determine if this is a nc or an ip sample */
        FILE** wigs = NULL;
        if ( NULL == strstr( argv[i+1], ".nc" ) )
        {    
            wigs = ip_wigs;
            ip_cntr += 1;
            cntr = ip_cntr;
        } else { 
            wigs = nc_wigs;
            nc_cntr += 1;
            cntr = nc_cntr;
        }
        
        wigs[cntr] = fopen( argv[i+1], "r" );
        if( NULL == wigs[cntr] )
        {
            perror( "FATAL     : Could not open wiggle file for reading" );
            exit( -1 );
        }
    }

    call_peaks_from_wiggles( ip_wigs, nc_wigs, (argc-1)/2, stdout, -0.9 );
    
    for( i = 0; i < (argc-1)/2; i++ )
    {
        fclose( ip_wigs[i] );
        fclose( nc_wigs[i] );
    }

    free( ip_wigs );
    free( nc_wigs );
    
    return 0;
}


