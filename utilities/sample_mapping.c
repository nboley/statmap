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
#include "../src/config_parsing.h"

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

    /* change to the output directory */
    int error = chdir( argv[1] );
    if( -1 == error )
    {
        perror( "FATAL       :  Cannot move into output directory ");
        exit( 1 );
    }
    
    /* read the assay info */
    struct args_t* args;
    read_config_file_from_disk( &args );
    
    /**** declare a few global variables that we will need ****/
    /* Seed the random number generator */
    srand ( time(NULL) );
    enum bool use_random_start = true;
    float max_prb_change_for_convergence = MAX_PRB_CHANGE_FOR_CONVERGENCE;
    if( args->num_threads == -1 )
    {
        /* try to get the number of available threads from the os */
        num_threads = get_nprocs();
        /* if we cant determine the number of threads, set it to 1 */
        if( num_threads <= 0 )
            num_threads = 1;
        
        /* never set the number of threads to more than 8, by default */
        if( num_threads > 8 )
            num_threads = 8;
        
        fprintf(stderr, "NOTICE      :  Number of threads is being set to %i \n", num_threads);
    } else {
        num_threads = args->num_threads;
    }
        
    min_num_hq_bps = args->min_num_hq_bps;
    /* END parse arguments */
    
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

    /* if necessary, load the fragment length distribution estimate */    
    if( args->assay_type == CHIP_SEQ ) {
        FILE* frag_len_fp = open_check_error( "estimated_fl_dist.txt", "r" );
        build_fl_dist_from_file( IP_mpd_rdb, frag_len_fp );
        fclose(frag_len_fp);
        
        build_chipseq_bs_density( IP_mpd_rdb->fl_dist );
        assert( IP_mpd_rdb->fl_dist != NULL );    
    }
    
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
        if( args->assay_type == CHIP_SEQ )
        {
            FILE* frag_len_fp = open_check_error( "estimated_fl_dist.txt", "r" );
            build_fl_dist_from_file( mpd_NC_rdb, frag_len_fp );
            fclose(frag_len_fp);
            
            build_chipseq_bs_density( mpd_NC_rdb->fl_dist );
        }
    }
    /*** END negative control data */    
    
    /* open the meta data file, to append to */
    FILE* meta_info_fp = open_check_error(RELAXED_SAMPLES_META_INFO_FNAME, "a");
    /* open the meta data file, to append to */
    FILE* ss_meta_info_fp = open_check_error(STARTING_SAMPLES_META_INFO_FNAME, "a");
    
    /** determine the next sample index **/
    int sample_index = determine_next_sample_index();

  
    /* actually take the samples */

    if( args->assay_type == CHIP_SEQ )
    {
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
        }
        else {
            take_chipseq_sample(
                IP_mpd_rdb, genome,
                sample_index, 
                ss_meta_info_fp, meta_info_fp,
                MAX_NUM_EM_ITERATIONS,
                max_prb_change_for_convergence
            );      
        }
    } else if ( args->assay_type == CAGE ){
        take_cage_sample(
            IP_mpd_rdb, genome,
            sample_index, 
            ss_meta_info_fp, meta_info_fp,
            MAX_NUM_EM_ITERATIONS,
            max_prb_change_for_convergence
        );
    } else {
        fprintf( stderr, "FATAL      : Unsupported assay '%i'\n", args->assay_type );
        exit( 1 );
    }    
    
    goto cleanup;

cleanup:
    fclose( meta_info_fp );
    
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
