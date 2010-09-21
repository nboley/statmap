#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>

/* reading directory entries */
#include <dirent.h>

/* to find out the  number of available processors */
#include <sys/sysinfo.h>

int num_threads = -1;
int min_num_hq_bps = -1;

#include "utility_common.h"

#include "../src/statmap.h"
#include "../src/trace.h"
#include "../src/wiggle.h"
#include "../src/sam.h"
#include "../src/fragment_length.h"
#include "../src/genome.h"
#include "../src/iterative_mapping.h"
#include "../src/mapped_read.h"
#include "../src/rawread.h"


struct fragment_length_dist_t* global_fl_dist;

void usage()
{
    fprintf( stderr, "Usage: ./build_wiggles_from_mapped_reads output_directory" );
    exit( 1 );
}

int main( int argc, char** argv )
{
    /* parse arguments */
    if( argc != 2  )
    {
        usage();
        exit(1);
    }

    /* Seed the random number generator */
    srand ( time(NULL) );
        
    enum assay_type_t assay_type = CHIP_SEQ;
    
    num_threads = 1;
    min_num_hq_bps = 10;
    /* END parse arguments */
        
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
    fprintf(stderr, "NOTICE      :  mmapping mapped reads DB.\n" );
    mmap_mapped_reads_db( mpd_rdb );
    fprintf(stderr, "NOTICE      :  indexing mapped reads DB.\n" );
    index_mapped_reads_db( mpd_rdb );
    
    /* load the fragment length distribution estimate */
    FILE* frag_len_fp = open_check_error( "estimated_fl_dist.txt", "r" );
    build_fl_dist_from_file( mpd_rdb, frag_len_fp );
    fclose(frag_len_fp);

    if( assay_type == CHIP_SEQ )
        build_chipseq_bs_density( mpd_rdb->fl_dist );
    
    /* load the rawreads */
    struct rawread_db_t* raw_rdb;
    populate_rawread_db( &raw_rdb, "reads.unpaired", "reads.pair1", "reads.pair2" );
    if( NULL == raw_rdb )
    {
        fprintf( stderr, "FATAL       :  Can not load raw reads.\n" );
        exit( 1 );
    }

    char* scaffold_names[2] = {"FWD", "BKWD"};
    
    /* initialize the traces to write to */
    struct trace_t* trace = NULL;
    init_trace( genome, &trace, 2, scaffold_names );
    
    bootstrap_traces_from_mapped_reads(
        mpd_rdb, trace, 
        update_chipseq_trace_expectation_from_location
    );    
    
    write_wiggle_from_trace_to_stream(
        trace, stdout, 0.01 );
    goto cleanup;

cleanup:
    free_genome( genome );

    close_mapped_reads_db( mpd_rdb );
    close_rawread_db( raw_rdb );
    
    return 0;
}
