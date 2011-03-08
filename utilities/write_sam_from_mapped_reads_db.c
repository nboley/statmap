#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "utility_common.h"

#include "../src/mapped_read.h"
#include "../src/genome.h"
#include "../src/sam.h"
#include "../src/rawread.h"
#include "../src/config_parsing.h"

/* make it possible to link */
int min_num_hq_bps = -1;
int num_threads = 1;

void usage()
{
    fprintf( stderr, "Usage: ./mapped_reads_to_sam output_directory [sample_number [use IP samples]]\n" );
}

int main( int argc, char** argv )
{
    if( argc != 2 && argc != 3 && argc != 4 )
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
    
    fprintf( stderr, "ASSAY TYPE  :  %i\n", args->assay_type );
    fprintf( stderr, "NC MPD RDB  :  %p\n", args->NC_rdb );

    /* default to using the marginal distribution, encoded with sample number 0 */
    int sample_num = -1;
    int USE_IP = -1;
    
    /* if this is an experiment with a negative control */
    if( NULL != args->NC_rdb )
    {
        /* if there is 1 additional paramater, assume we want the marginal */
        if( argc == 3 ) {
            sample_num = 0;
            USE_IP = atoi( argv[2] );
        } else 
        /* we've specified a sample number */
        {
            sample_num = atoi( argv[2] );
            USE_IP = atoi( argv[3] );
        }

        if( USE_IP < 0 ) {
            usage();
            fprintf( stderr, "ERROR    : Invalid use IP flag (%i) - expecting 0 ( for use negative control ) or 1 ( for use ip )\n", sample_num );
            fprintf( stderr, "HINT     : This experiment used a negative control - you need to specify which reads to output.\n", sample_num );
            exit(1);
        }
    } else {
        if( argc == 3 ) {
            sample_num = atoi( argv[2] );
        }
    }
    
    if( sample_num < 0  ) {
        usage();
        fprintf( stderr, "ERROR    : Invalid Sample Number %i\n", sample_num );
        exit(1);
    }
        
    /* Load the genome */
    struct genome_data* genome;
    load_genome_from_disk( &genome, GENOME_FNAME );

    /* load the mapped read db */
    char* mpd_rd_fname = NULL;
    if( USE_IP == 0 ) {
        mpd_rd_fname = "mapped_NC_reads.db";
    } else {
        mpd_rd_fname = "mapped_reads.db";
    }
    struct mapped_reads_db* mpd_rdb;
    open_mapped_reads_db( &mpd_rdb, mpd_rd_fname );

    /* load the rawreads */
    struct rawread_db_t* raw_rdb;
    if( USE_IP == 0 ) {
        populate_rawread_db( &raw_rdb, "reads.NC.unpaired", "reads.NC.pair1", "reads.NC.pair2" );
    } else {
        populate_rawread_db( &raw_rdb, "reads.unpaired", "reads.pair1", "reads.pair2" );
    }
    
    if( NULL == raw_rdb )
    {
        fprintf( stderr, "FATAL       :  Can not load raw reads.\n" );
        exit( 1 );
    }    

    /* default to reset conditional probabilities */
    enum bool reset_cond_prbs = true;
    
    /* If we specified a marginal desnity */
    if( sample_num > 0 ) {
        /* tell the sam printing code to not reset the cond prbs */
        reset_cond_prbs = false;
        
        /* first, load the fl dist */
        FILE* fl_dist_fp = NULL;
        fl_dist_fp = fopen( "estimated_fl_dist.txt", "r"  );
        init_fl_dist_from_file( &(mpd_rdb->fl_dist), fl_dist_fp );
        fclose( fl_dist_fp  );
        build_chipseq_bs_density( mpd_rdb->fl_dist );
        
        /* load the trace that stores the marginal read 
           density that we are interested in */
        /* this is a bit confusing because with assay/experiment combos
           that use a negative control, there are 2 mappings for each experiment.
           So, we check the command line paramater specifying the negative control.
           If it is set, we assuem the experiment had a negative control and load
           it depending on the flag's value.
        */
        struct trace_t* traces;
        char traces_fname[500];
        if( -1 == USE_IP ) {
            sprintf( traces_fname, "./samples/sample%i.bin.trace", sample_num );
        } else if( 1 == USE_IP ) {
            sprintf( traces_fname, "./samples/sample%i.ip.bin.trace", sample_num );
        } else if( 0 == USE_IP ) {
            sprintf( traces_fname, "./samples/sample%i.nc.bin.trace", sample_num );
        }
        
        fprintf( stderr, "TRACE FNAME :  %s\n", traces_fname );
        load_trace_from_file( &traces, traces_fname  );

        /* update the read conditional probabilities based upon the assay
           and the correct trace */
        mmap_mapped_reads_db( mpd_rdb );
        index_mapped_reads_db( mpd_rdb );
        update_cond_prbs_from_trace_and_assay_type( 
            mpd_rdb, traces, genome, args->assay_type );

        close_traces( traces );
    }
        
    write_mapped_reads_to_sam( 
        raw_rdb, mpd_rdb, genome, reset_cond_prbs, false, stdout );
    
    goto cleanup;
    
cleanup:
    free( args );
    free_genome( genome );
    close_mapped_reads_db( mpd_rdb );
    close_rawread_db( raw_rdb );    

    return 0;
}



