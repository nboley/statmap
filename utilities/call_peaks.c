#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int num_threads = -1;
int min_num_hq_bps = -1;

#include "../src/wiggle.h"
#include "../src/trace.h"
#include "../src/genome.h"

void usage()
{
    fprintf(stderr, "Usage: ./call_peaks genome.bin IP_bs_sample.wig NC_bs_sample.wig > called_peaks.wig\n");
}


int 
main( int argc, char** argv )
{
    if( argc != 4 )
    {
        usage();
        exit( -1 );
    }

    /* Load the genome */
    struct genome_data* genome;
    load_genome_from_disk( &genome, argv[1] );
    
    /* open the wiggles */
    FILE* ip_wig = fopen( argv[2], "r" );
    if( NULL == ip_wig )
    {
        perror( "FATAL     : Could not open IP wiggle file for reading" );
        exit( -1 );
    }
    
    FILE* nc_wig = fopen( argv[3], "r" );
    if( NULL == nc_wig )
    {
        perror( "FATAL     : Could not open NC wiggle file for reading" );
        exit( -1 );
    }
    
    /* init the trace */
    struct trace_t* trace;
    /* two tracks correspond to pos and neg strands */
    init_traces( genome, &trace, 2 );
    
    /* update the trace */
    call_peaks_from_wiggles( ip_wig, nc_wig, genome, trace );
    
    fclose( ip_wig );
    fclose( nc_wig );
    
    return 0;
}


