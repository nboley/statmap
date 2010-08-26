#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>

/* reading directory entries */
#include <dirent.h>

int num_threads = -1;
int min_num_hq_bps = -1;

#include "../src/config.h"
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

int
determine_next_sample_index()
{
    DIR *dp;
    struct dirent *ep;

    dp = opendir ( BOOTSTRAP_SAMPLES_ALL_PATH );
    if( NULL == dp )
    {
        perror( "FATAL       : Could not open bootstrap samples directory." );
        exit( 1 );
    }

    int sample_index = 0;
    while( 0 != ( ep = readdir (dp) ) )
    {
        /* skip the references */
        if( ep->d_name[0] == '.' )
            continue;     
        sample_index += 1;
    }

    return sample_index;
}

void usage()
{
    fprintf( stderr, "Usage: ./sample_mapping.c output_directory\n" );
}

void
populate_rawread_db( 
    struct rawread_db_t** raw_rdb,
    char* unpaired_fname,
    char* pair1_fname,
    char* pair2_fname )
{
    *raw_rdb = NULL;
    
    /* try to open the single end file to see if the reads are single ended */
    FILE* tmp = fopen( unpaired_fname, "r" );
    if( tmp != NULL )
    {
        fclose( tmp );
        init_rawread_db( raw_rdb );
        add_single_end_reads_to_rawread_db(
            *raw_rdb, unpaired_fname, FASTQ 
        );
        
        return;
    } 

    tmp = fopen( pair1_fname, "r" );
    if( tmp != NULL )
    {
        fclose( tmp );
        init_rawread_db( raw_rdb );
        
        add_paired_end_reads_to_rawread_db(
            *raw_rdb, 
            pair1_fname,
            pair2_fname,
            FASTQ 
        );
        
        return;
    }
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

    close_mapped_reads_db( IP_mpd_rdb );
    close_rawread_db( raw_rdb );

    if( NULL != raw_NC_rdb )
    {
        close_rawread_db( raw_NC_rdb );
        close_mapped_reads_db( mpd_NC_rdb );
    }
    
    return 0;
}
