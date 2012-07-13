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
    fprintf( stderr, "Usage: ./mapped_reads_to_sam output_directory input.bin.trace\n" );
}

int main( int argc, char** argv )
{
    if( argc != 3 )
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

    /* load the mapped read db */
    char* mpd_rd_fname = NULL;
    mpd_rd_fname = "mapped_reads.db";
    struct mapped_reads_db* mpd_rdb;
    open_mapped_reads_db( &mpd_rdb, mpd_rd_fname );

    struct cond_prbs_db_t* cond_prbs_db;
    init_cond_prbs_db_from_mpd_rdb( &cond_prbs_db, mpd_rdb );
    
    /* load the rawreads */
    struct rawread_db_t* raw_rdb;
    populate_rawread_db( &raw_rdb, "reads.unpaired", "reads.pair1", "reads.pair2" );
    
    if( NULL == raw_rdb )
    {
        fprintf( stderr, "FATAL       :  Can not load raw reads.\n" );
        exit( 1 );
    }    

    /* default not to reset conditional probabilities */
    enum bool reset_cond_prbs = false;
    
    /* load the trace that stores the marginal read 
       density that we are interested in */
    struct trace_t* traces;
    load_trace_from_file( &traces, argv[2] );

    /* update the read conditional probabilities based upon the assay
       and the correct trace */
    mmap_mapped_reads_db( mpd_rdb );
    index_mapped_reads_db( mpd_rdb );
    update_cond_prbs_from_trace_and_assay_type( 
        mpd_rdb, cond_prbs_db, traces, genome, args->assay_type );

    close_traces( traces );
        
    write_mapped_reads_to_sam( 
        raw_rdb, mpd_rdb, cond_prbs_db, genome, reset_cond_prbs, false, stdout );
    
    goto cleanup;
    
cleanup:
    free( args );
    free_genome( genome );
    close_mapped_reads_db( &mpd_rdb );
    close_rawread_db( raw_rdb );    

    return 0;
}


