#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>

int num_threads = -1;
int min_num_hq_bps = -1;

#include "../src/wiggle.h"
#include "../src/trace.h"
#include "../src/genome.h"

static FILE* 
open_check_error( char* fname, char* file_mode )
{
    FILE* tmp;
    tmp = fopen( fname, file_mode );
    if( tmp == NULL )
    {
        fprintf( stderr, "Error opening '%s\n'", fname );
        exit( -1 );
    }
    return tmp;
}

void usage()
{
    fprintf(stderr, "Usage: ./call_peaks output_directory > called_peaks.wig\n");
}


int 
main( int argc, char** argv )
{
    if( argc != 2 )
    {
        usage();
        exit( -1 );
    }

    char* curr_wd = NULL;
    curr_wd = getcwd( NULL, 0 );
    if( curr_wd == NULL )
    {
        perror( "FATAL       :  Cannot get current working directory ");
        exit( 1 );        
    }
    
    /* change to the output directory */
    int error = chdir( argv[1] );
    if( -1 == error )
    {
        perror( "FATAL       :  Cannot move into output directory ");
        exit( 1 );
    }
    
    /* Load the genome */
    struct genome_data* genome;
    load_genome_from_disk( &genome, GENOME_FNAME );
    
    /* init the trace */
    struct trace_t* trace;
    /* two tracks correspond to pos and neg strands */
    init_traces( genome, &trace, 2 );
    
    int si = 0; /* sample index */
    int bi; /* bootstrap index */
    for( bi = 0; bi < NUM_BOOTSTRAP_SAMPLES; bi++ )
    {
        char fname[500];

        /* open the IP wiggle file */
        sprintf( fname, "%s/sample%i/bssample%i.ip.wig", 
                 BOOTSTRAP_SAMPLES_ALL_PATH, si+1, bi+1 );
        FILE* ip_wig = open_check_error( fname, "r" ); 

        /* open the NC wiggle file */
        sprintf( fname, "%s/sample%i/bssample%i.nc.wig", 
                 BOOTSTRAP_SAMPLES_ALL_PATH, si+1, bi+1 );
        FILE* nc_wig = open_check_error( fname, "r" ); 

        call_peaks_from_wiggles( ip_wig, nc_wig, genome, trace );

        fclose( ip_wig );
        fclose( nc_wig );
    }
    
    char* track_names[2] = {"fwd_strand_peaks", "bkwd_strnad_peaks"};
    write_wiggle_from_trace( trace, genome->chr_names, track_names, "test_out.wig", 0.0 );
    
    close_traces( trace );
    free_genome( genome );
    free( curr_wd );
    
    return 0;
}


