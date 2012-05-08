#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <float.h>
#include <assert.h>

#include "../src/config_parsing.h"
#include "../src/genome.h"
#include "../src/mapped_read.h"
#include "../src/iterative_mapping.h"
#include "../src/trace.h"
#include "../src/trace_agg_fns.h"

/* make it possible to link */
int min_num_hq_bps = -1;
int num_threads = 2;

void usage()
{
    fprintf(stderr, "Usage: ./build_min_trace (min|max|sum) output.bin.trace statmap_output_directory/ sample_number\n");
}

void
bootstrap_trace(
    struct trace_t*         trace,
    struct mapped_reads_db* mpd_rdb,
    struct cond_prbs_db_t*  cond_prbs_db,
    void (* const update_trace_expectation_from_location)(
        const struct trace_t* const traces, 
        const struct mapped_read_location* const loc,
        float cond_prob
    )
)
{
    bootstrap_traces_from_mapped_reads(
        mpd_rdb, cond_prbs_db, trace, update_trace_expectation_from_location
    );
}

int
main(int argc, char** argv)
{
    if( argc != 5 )
    {
        usage();
        exit( -1 );
    }

    /*** determine the aggregate type ***/
    TRACE_TYPE (*agg_fn)( const TRACE_TYPE, const TRACE_TYPE ) = NULL;

    if( 0 == strcmp("max", argv[1] ) )
    {
        agg_fn = max;
    } else if( 0 == strcmp("min", argv[1] ) ) {
        agg_fn = min;
    } else if( 0 == strcmp("sum", argv[1] ) ) {
        agg_fn = sum;
    } else {
        fprintf( stderr, "FATAL     : Unrecognized aggregate '%s'\n", argv[1] );
        exit( -1 );
    }

    /*** load data structures and metadata from statmap output directory ***/

    // change working directory to the statmap output directory
    int rv;
    rv = chdir( argv[3] );
    if( rv != 0 )
    {
        perror( "FATAL : Cannot move into output directory" );
        exit( -1 );
    }

    // load saved args from config.dat
    struct args_t* args;
    read_config_file_fname_from_disk( "config.dat", &args );

    int sample_num = atoi( argv[4] );
    int USE_IP = -1;

    /* determine assay type and # of starting locations from configuration file
     * (so we can set update_trace_expectation_from_location accordingly) */
    enum assay_type_t assay_type = args->assay_type;
    assert( assay_type == CAGE || assay_type == CHIP_SEQ );

    /* set the update_trace_expectation_from_location function ptr
     * (TODO: and the track names) based on the assay type */
    void (* update_trace_expectation_from_location)(
        const struct trace_t* const traces, 
        const struct mapped_read_location* const loc,
        float cond_prob
    ) = NULL;
    assert( assay_type == CAGE || assay_type == CHIP_SEQ );
    if( assay_type == CAGE )
    {
        update_trace_expectation_from_location = \
            &update_CAGE_trace_expectation_from_location;
    } else if( assay_type == CHIP_SEQ )
    {
        update_trace_expectation_from_location = \
            &update_chipseq_trace_expectation_from_location;
    }

    /* define possible trace track names (choose 1 depending on assay) */
    // TODO - handle all possibilities for chipseq traces
    char* track_names[2] = { "fwd_strnd_read_density", "rev_strnd_read_density" };
    char* ip_track_names[2] = {"IP_fwd_strand", "IP_bkwd_strand"};
    char* nc_track_names[2] = {"NC_fwd_strand", "NC_bkwd_strand"};

    /* load the genome */
    struct genome_data* genome;
    load_genome_from_disk( &genome, GENOME_FNAME );

    /* load the mapped reads db */
    struct mapped_reads_db* mpd_rdb = NULL;
    open_mapped_reads_db( &mpd_rdb, MAPPED_READS_DB_FNAME ); // XXX: different for chipseq?
    mmap_mapped_reads_db( mpd_rdb );
    index_mapped_reads_db( mpd_rdb );

    /* initialize the cond prbs db */
    struct cond_prbs_db_t* cond_prbs_db = NULL;
    init_cond_prbs_db_from_mpd_rdb( &cond_prbs_db, mpd_rdb );

    /*** bootstrap traces ***/

    /* build the first trace as the trace to update - if sample_number was specified,
     * build from those cond probs*/
    struct trace_t* update_trace;
    //init_trace( genome, &update_trace, 2, track_names );

    /* load the trace that stores the marginal read 
       density that we are interested in */
    /* this is a bit confusing because with assay/experiment combos
       that use a negative control, there are 2 mappings for each experiment.
       So, we check the command line paramater specifying the negative control.
       If it is set, we assuem the experiment had a negative control and load
       it depending on the flag's value.
    */
    char update_trace_fname[500];
    if( -1 == USE_IP ) {
        sprintf( update_trace_fname, "./samples/sample%i.bin.trace", sample_num );
    } else if( 1 == USE_IP ) {
        sprintf( update_trace_fname, "./samples/sample%i.ip.bin.trace", sample_num );
    } else if( 0 == USE_IP ) {
        sprintf( update_trace_fname, "./samples/sample%i.nc.bin.trace", sample_num );
    }
        
    load_trace_from_file( &update_trace, update_trace_fname  );

    // update cond prbs from specified trace
    update_cond_prbs_from_trace_and_assay_type( mpd_rdb,
                                                cond_prbs_db,
                                                update_trace,
                                                genome,
                                                assay_type );

    fprintf( stderr, "Bootstrapping %i samples.", NUM_BOOTSTRAP_SAMPLES );

    /* build min aggregate, two traces at a time */
    int i;
    // TODO - NUM_BOOTSTRAP_SAMPLES as parameter
    struct trace_t* curr_trace;
    init_trace( genome, &curr_trace, 2, track_names );
    for( i = 1; i < NUM_BOOTSTRAP_SAMPLES; i++ )
    {
        if(  i%(NUM_BOOTSTRAP_SAMPLES/10)  == 0 )
            fprintf( stderr, " %.1f%%...", (100.0*i)/NUM_BOOTSTRAP_SAMPLES );

        /* bootstrap current trace */
        /* note - bootstrap_trace calls bootstrap_traces_from_mapped_reads,
         * which zeros out curr_trace */
        bootstrap_trace( curr_trace,
                         mpd_rdb,
                         cond_prbs_db,
                         update_trace_expectation_from_location 
                       );

        /* aggregate curr_trace into update_trace with agg_fn */
        aggregate_over_trace_pairs( update_trace, curr_trace, agg_fn );
    }
    close_traces( curr_trace );

    write_trace_to_file( update_trace, argv[2] );

    close_traces( update_trace );

    return 0;
}
