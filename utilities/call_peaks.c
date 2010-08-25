#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "../src/config.h"
#include "../src/genome.h"
#include "../src/peak_calling.h"

void usage()
{
    fprintf(stderr, "Usage: ./call_peaks output_directory > called_peaks.wig\n");
}

int num_threads = -1;
int min_num_hq_bps = -1;

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

    call_peaks( genome );

    free_genome( genome );
    free( curr_wd );
    
    return 0;
}


