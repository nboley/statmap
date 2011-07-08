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

    struct mapped_reads_db* new_mpd_rdb;
    new_mapped_reads_db( &new_mpd_rdb, "new_mapped_reads.db" );

    struct mapped_read_t* mapped_rd;
    
    error = get_next_read_from_mapped_reads_db( 
        mpd_rdb, &mapped_rd
    );

    while( NULL != mapped_rd )
    {
        // printf("%i\t%i\n", mapped_rd->read_id, mapped_rd->num_mappings );
        unsigned int i;
        enum bool discovered_error = false;
        for( i = 0; i < mapped_rd->num_mappings; i++ )
        {
            // Make sure that the read is actually in the genome
            int chr = get_chr_from_mapped_read_location( mapped_rd->locations + i );
            int start = get_start_from_mapped_read_location( mapped_rd->locations + i );
            int stop = get_stop_from_mapped_read_location( mapped_rd->locations + i );
            if( chr > genome->num_chrs
                || start < 0 || stop < 0
                || (unsigned int) stop >= genome->chr_lens[ chr ]
                )
            {
                printf("\t\t%i (%i)\t%i-%i ( %i )\n", chr, genome->num_chrs, start, stop, genome->chr_lens[ chr ] );
                discovered_error = true;
            } 
        }

        /* if there were no errors, add the read to the new rdb */
        if( discovered_error == false )
        {
            add_read_to_mapped_reads_db( new_mpd_rdb, mapped_rd );
        }
        
        // print_mapped_locations( mapped_rd->locations );
        
        free_mapped_read( mapped_rd );
                
        error = get_next_read_from_mapped_reads_db( 
            mpd_rdb, &mapped_rd
        );

    }
    
    close_mapped_reads_db( &new_mpd_rdb );
    close_mapped_reads_db( &mpd_rdb );
    
    return 0;
}
