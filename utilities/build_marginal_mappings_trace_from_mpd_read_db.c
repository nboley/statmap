#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include <unistd.h>

#include "../src/config.h"
#include "../src/genome.h"
#include "../src/trace.h"
#include "../src/iterative_mapping.h"
#include "../src/mapped_read.h"

/* make it possible to link */
int min_num_hq_bps = -1;
int num_threads = 1;

/* forward declarations */
struct mapped_reads_db_t;

void usage()
{
    fprintf(stderr, "Usage: ./build_marginal_mappings_trace_from_mpd_read_db output_directory [filter_threshold] > output.wig\n");
}


int main( int argc, char** argv )
{
    /* parse the arguments */
    if( argc <= 1 || argc >= 4 )
    {
        usage();
        exit( -1 );
    }

    double filter_threshold = 0;
    if( argc == 3 )
    {
        filter_threshold = atof( argv[2] );
    }

    /* change into the output directory */
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
    
    /* load the mapped reads */
    char* mpd_rd_fname = "mapped_reads.db";
    struct mapped_reads_db* mpd_rdb;
    open_mapped_reads_db( &mpd_rdb, mpd_rd_fname );    
    struct cond_prbs_db_t* cond_prbs_db;
    init_cond_prbs_db_from_mpd_rdb( &cond_prbs_db, mpd_rdb );
    reset_all_read_cond_probs( mpd_rdb, cond_prbs_db );

    /* initialize a trace */
    struct trace_t* trace;
    char* track_names[2] = {"fwd_strand", "bkwd_strand"};
    init_trace( genome, &trace, 2, track_names );
    
    /*** update the trace from the naive kernel ***/
    /* we hack this a bit by noting that the naive kernel is just the 
       trace kernel ( stranded, which we nearly always want ) */
    update_traces_from_mapped_reads( 
        mpd_rdb, cond_prbs_db, trace, update_CAGE_trace_expectation_from_location);
    
    /* write the trace to stdout */
    write_wiggle_from_trace_to_stream( trace, stdout, filter_threshold );

    /* cleanup */
    goto cleanup;
    
cleanup:
    free_genome( genome );
    close_mapped_reads_db( &mpd_rdb );
    
}
