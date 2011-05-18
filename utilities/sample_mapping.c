#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>

/* to find out the  number of available processors */
#include <sys/sysinfo.h>

int num_threads = -1;
int min_num_hq_bps = -1;

#include "utility_common.h"

#include "../src/statmap.h"
#include "../src/config_parsing.h"
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
    fprintf( stderr, "Usage: ./sample_mapping.c output_directory\n" );
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
    enum bool use_random_start = true;
    float max_prb_change_for_convergence = 1e-2;

    /* try to get the number of available threads from the os */
    num_threads = MAX( 8, get_nprocs() );
    /* if we cant determine the number of threads, set it to 1 */
    if( num_threads <= 0 )
        num_threads = 1;
    fprintf(stderr, "NOTICE      :  Number of threads is being set to %i \n", num_threads);

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
    struct mapped_reads_db* IP_mpd_rdb;
    open_mapped_reads_db( &IP_mpd_rdb, mpd_rd_fname );
    fprintf(stderr, "NOTICE      :  mmapping mapped reads DB.\n" );
    mmap_mapped_reads_db( IP_mpd_rdb );
    fprintf(stderr, "NOTICE      :  indexing mapped reads DB.\n" );
    index_mapped_reads_db( IP_mpd_rdb );
    
    /* load the fragment length distribution estimate */
    FILE* frag_len_fp = open_check_error( "estimated_fl_dist.txt", "r" );
    build_fl_dist_from_file( IP_mpd_rdb, frag_len_fp );
    fclose(frag_len_fp);

    if( assay_type == CHIP_SEQ )
        build_chipseq_bs_density( IP_mpd_rdb->fl_dist );

    assert( IP_mpd_rdb->fl_dist != NULL );
    
    /* load the rawreads */
    struct rawread_db_t* raw_rdb;
    populate_rawread_db( &raw_rdb, "reads.unpaired", "reads.pair1", "reads.pair2" );
    if( NULL == raw_rdb )
    {
        fprintf( stderr, "FATAL       :  Can not load raw reads.\n" );
        exit( 1 );
    }

    /**** 
     **** Try loading the negative control data. If we can load the raw 
     **** data, then assume the mapped reads exist.
     ****
     ****/
    struct rawread_db_t* raw_NC_rdb = NULL;
    populate_rawread_db( 
        &raw_NC_rdb, "reads.NC.unpaired", "reads.NC.pair1", "reads.NC.pair2" );
    /* we dont check for NULL, because there may not be negative control reads */
    
    /* if there are negative control reads,
        then assume there are mapped nc reads */
    struct mapped_reads_db* mpd_NC_rdb = NULL;
    if( NULL != raw_NC_rdb )
    {
        char* mpd_NC_rd_fname = "mapped_NC_reads.db";
        open_mapped_reads_db( &mpd_NC_rdb, mpd_NC_rd_fname );
        fprintf(stderr, "NOTICE      :  mmapping mapped reads DB.\n" );
        mmap_mapped_reads_db( mpd_NC_rdb );
        fprintf(stderr, "NOTICE      :  indexing mapped reads DB.\n" );
        index_mapped_reads_db( mpd_NC_rdb );

        /* load the fragment length distribution estimate */
        FILE* frag_len_fp = open_check_error( "estimated_fl_dist.txt", "r" );
        build_fl_dist_from_file( mpd_NC_rdb, frag_len_fp );
        fclose(frag_len_fp);
        
        if( assay_type == CHIP_SEQ )
            build_chipseq_bs_density( mpd_NC_rdb->fl_dist );

    }
    
    
    /* open the meta data file, to append to */
    FILE* meta_info_fp = open_check_error(RELAXED_SAMPLES_META_INFO_FNAME, "a");
    
    /** determine the next sample index **/
    int sample_index = determine_next_sample_index();
    
    /*** END negative control data */
    if( NULL != mpd_NC_rdb )
    {
        take_chipseq_sample_wnc(
            IP_mpd_rdb, mpd_NC_rdb,
            genome,
            meta_info_fp,
            sample_index,
            max_prb_change_for_convergence,
            use_random_start 
        );
    } else {
        /*
        error = sample_relaxed_mapping(  
            raw_rdb, mpd_rdb, genome,
            assay_type, use_random_start, max_prb_change_for_convergence, 
            argv[2], sam_ofname
        );
        */
    }
    
    fclose( meta_info_fp );
    
    
    
    goto cleanup;

cleanup:
    free_genome( genome );

    close_mapped_reads_db( &IP_mpd_rdb );
    close_rawread_db( raw_rdb );

    if( NULL != raw_NC_rdb )
    {
        close_rawread_db( raw_NC_rdb );
        close_mapped_reads_db( &mpd_NC_rdb );
    }
    
    return 0;
}
