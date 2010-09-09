#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "utility_common.h"

#include "../src/mapped_read.h"
#include "../src/genome.h"
#include "../src/sam.h"
#include "../src/rawread.h"

/* make it possible to link */
int min_num_hq_bps = -1;

void usage()
{
    fprintf( stderr, "Usage: ./mapped_reads_to_sam output_directory\n" );
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

    /* Load the genome */
    struct genome_data* genome;
    load_genome_from_disk( &genome, GENOME_FNAME );
    
    /* load the mapped read db */
    char* mpd_rd_fname = "mapped_reads.db";
    struct mapped_reads_db* mpd_rdb;
    open_mapped_reads_db( &mpd_rdb, mpd_rd_fname );

    /* load the rawreads */
    struct rawread_db_t* raw_rdb;
    populate_rawread_db( &raw_rdb, "reads.unpaired", "reads.pair1", "reads.pair2" );
    if( NULL == raw_rdb )
    {
        fprintf( stderr, "FATAL       :  Can not load raw reads.\n" );
        exit( 1 );
    }
    
    write_mapped_reads_to_sam( raw_rdb, mpd_rdb, genome, true, false, stdout );
    
    goto cleanup;
    
cleanup:
    free_genome( genome );
    close_mapped_reads_db( mpd_rdb );
    close_rawread_db( raw_rdb );

    return 0;
}



