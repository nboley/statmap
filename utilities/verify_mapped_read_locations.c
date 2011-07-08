#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "utility_common.h"

#include "../src/mapped_read.h"
#include "../src/genome.h"
#include "../src/sam.h"
#include "../src/rawread.h"
#include "../src/config_parsing.h"
#include "../src/fragment_length.h"
#include "../src/trace.h"
#include "../src/iterative_mapping.h"

/* make it possible to link */
int min_num_hq_bps = -1;
int num_threads = 1;

void usage()
{
    fprintf( stderr, "Usage: ./verify_mapped_read_locations output_directory\n" );
}

int main( int argc, char** argv )
{
    if( argc != 2 )
    {
        usage();
        exit(1);
    }

    /* change to the output directory */
    int error = chdir( argv[1] );
    if( -1 == error )
    {
        perror( "FATAL       :  Cannot move into output directory ");
        exit( 1 );
    }
    
    /* load the configuration data */
    struct args_t* args;
    read_config_file_from_disk( &args );

    /* Load the genome */
    struct genome_data* genome;
    load_genome_from_disk( &genome, GENOME_FNAME );

    
    char* mpd_rd_fname = "mapped_reads.db";
    struct mapped_reads_db* mpd_rdb;
    open_mapped_reads_db( &mpd_rdb, mpd_rd_fname );
    rewind_mapped_reads_db( mpd_rdb );

    struct mapped_read_t* mapped_rd;
    
    error = get_next_read_from_mapped_reads_db( 
        mpd_rdb, &mapped_rd
    );

    while( NULL != mapped_rd )
    {
        printf("%i\n", mapped_rd->read_id );
        
        free_mapped_read( mapped_rd );
                
        error = get_next_read_from_mapped_reads_db( 
            mpd_rdb, &mapped_rd
        );

    }
    
    return 0;
}
