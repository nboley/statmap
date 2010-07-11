#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>

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

/* This is to silence a link error with including genome.o. It's OK because this never
   indexes the genome, so it never needs to free it */
void free_tree( void )
{ return; }

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
sample_relaxed_mapping(  
    struct rawread_db_t* rawread_db,
    struct mapped_reads_db* rdb, 
    struct genome_data* genome,

    enum assay_type_t assay_type,
    enum bool use_random_start,
    float max_prb_change_for_convergence,

    char* wiggle_ofname,
    char* sam_ofname
)
{
    assert( UNKNOWN != assay_type );
    
    clock_t start, stop;
    
    int error = 0;
    
    /* store the assay specific callbacks */
    void (*update_expectation)(
        const struct trace_t* const traces, 
        const struct mapped_read_location* const loc) 
        = NULL;
    
    struct update_mapped_read_rv_t 
        (*update_reads)( const struct trace_t* const traces, 
                         const struct mapped_read_t* const r  )
        = NULL;

    char** track_names = NULL;

    int trace_size = -1;

    switch( assay_type )
    {
    case CAGE:
        update_expectation = update_CAGE_trace_expectation_from_location;
        update_reads = update_CAGE_mapped_read_prbs;
        trace_size = 2;
        track_names = malloc( trace_size*sizeof(char*) );
        track_names[0] = "fwd_strnd_read_density"; 
        track_names[1] = "rev_strnd_read_density";
        break;
    
    case CHIP_SEQ:
        update_expectation = update_chipseq_trace_expectation_from_location;
        update_reads = update_chipseq_mapped_read_prbs;
        trace_size = 1;
        track_names = malloc( trace_size*sizeof(char*) );
        track_names[0] = "read_density";
        break;
    
    default:
        fprintf( stderr, "WARNING     :  Can not iteratively map for assay type '%u'.\n", assay_type);
        exit( 1 );
    }
    
    /* BUG!!! Set the global fl dist */
    global_fl_dist = rdb->fl_dist;
    
    /* reset the read cond prbs under a uniform prior */
    /* this should also cache the mapped reads databse */
    reset_all_read_cond_probs( rdb );
    
    struct trace_t* starting_trace;
    
    /* iteratively map from a uniform prior */
    start = clock();
    fprintf(stderr, "NOTICE      :  Starting iterative mapping.\n" );
    
    /* initialize the trace that we will store the expectation in */
    init_traces( genome, &starting_trace, trace_size );
    if( use_random_start == true )
    {
        fprintf(stderr, "DEBUG       :  Building Random Starting Trace.\n" );
        build_random_starting_trace( 
            starting_trace, rdb, 
            update_expectation,
            update_reads
        );
    } else {
        fprintf(stderr, "DEBUG       :  Building Uniform Starting Trace.\n" );
        set_trace_to_uniform( starting_trace, 1 );
    }
    
    fprintf(stderr, "DEBUG       :  Updating mapping.\n" );
    error = update_mapping (
        rdb, 
        starting_trace,
        MAX_NUM_EM_ITERATIONS,
        max_prb_change_for_convergence,
        update_expectation,
        update_reads
    );
    
    stop = clock();
    fprintf(stderr, "PERFORMANCE :  Maximized LHD in %.2lf seconds\n", 
            ((float)(stop-start))/CLOCKS_PER_SEC );

    /* write the mapped read density to a wiggle */
    if( NULL != wiggle_ofname )
    {
        start = clock();
        fprintf(stderr, "NOTICE      :  Writing relaxed mapped read density to wiggle file.\n" );
        
        write_wiggle_from_trace( 
            starting_trace, 
            genome->chr_names, track_names,
            wiggle_ofname, max_prb_change_for_convergence );
        
        stop = clock();
        fprintf(stderr, 
                "PERFORMANCE :  Wrote relaxed mapped read density to wiggle in %.2lf seconds\n", 
                ((float)(stop-start))/CLOCKS_PER_SEC );
    }
    
    /* write the mapped reads to SAM */
    if( NULL != sam_ofname )
    {
        start = clock();
        fprintf(stderr, "NOTICE      :  Writing relaxed mapped reads to SAM file.\n" );
    
        FILE* sam_ofp = fopen( sam_ofname, "w+" );
        write_mapped_reads_to_sam( 
            rawread_db, rdb, genome, false, false, sam_ofp );
        fclose( sam_ofp );    
        
        stop = clock();
        fprintf(stderr, 
                "PERFORMANCE :  Wrote relaxed mapped reads to sam in %.2lf seconds\n", 
                ((float)(stop-start))/CLOCKS_PER_SEC );
    }
    
    goto cleanup;

cleanup:

    free( track_names );

    close_traces( starting_trace );
    
    return 0;    



}


void usage()
{
    fprintf( stderr, "Usage: ./mapped_reads_to_sam output_directory output.wig [ output.sam ]\n" );
}

int main( int argc, char** argv )
{
    /* parse arguments */
    if( argc < 3 || argc > 4 )
    {
        usage();
        exit(1);
    }

    enum assay_type_t assay_type = CHIP_SEQ;
    enum bool use_random_start = false;
    float max_prb_change_for_convergence = 1e-2;
    num_threads = 1;
    min_num_hq_bps = 10;
    
    /* END parse arguments */

    char* curr_wd = NULL;
    curr_wd = getcwd( NULL, 0 );
    if( curr_wd == NULL )
    {
        perror( "FATAL       :  Cannot get current working directory ");
        exit( 1 );        
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
    init_genome( &genome );
    FILE* genome_fp = open_check_error( "genome.fa", "r" );
    add_chrs_from_fasta_file( genome,  genome_fp );

    // BUG Why does it segfault when I close the file handle?
    
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
    init_rawread_db( &raw_rdb );
    
    /* try to open the single end file to see if the reads are single ended */
    FILE* tmp = fopen( "reads.unpaired", "r" );
    if( tmp != NULL )
    {
        fclose( tmp );

        add_single_end_reads_to_rawread_db(
            raw_rdb, "reads.unpaired", FASTQ 
        );
    } 
    /* else, assume the reads are paired */
    else {
        add_paired_end_reads_to_rawread_db(
            raw_rdb, 
            "reads.pair1",
            "reads.pair2",
            FASTQ 
        );
    }

    /* change to the curr working directory */
    error = chdir( curr_wd );
    if( -1 == error )
    {
        perror( "FATAL       :  Cannot move back into the initial working directory ");
        exit( 1 );
    }
    free( curr_wd );

    char* sam_ofname = NULL;
    if( argc == 4 )
        sam_ofname = argv[3];

    error = sample_relaxed_mapping(  
        raw_rdb, mpd_rdb, genome,
        assay_type, use_random_start, max_prb_change_for_convergence, 
        argv[2], sam_ofname
    );
    
    goto cleanup;

cleanup:
    free_genome( genome );
    close_mapped_reads_db( mpd_rdb );
    close_rawread_db( raw_rdb );

    return 0;
}
