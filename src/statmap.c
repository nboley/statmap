/* Copyright (c) 2009-2012, Nathan Boley */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h> /* gettimeofday() */
#include <libgen.h> /* get the base directory name */

#include <string.h>
#include <limits.h>

#include <pthread.h>

/* to find out the  number of available processors */
#include <sys/sysinfo.h>
/* mkdir */
#include <sys/stat.h>
/* chdir */
#include <sys/unistd.h>
/* check for system signals */
#include <signal.h>
#include <sys/wait.h>

#include "config_parsing.h"
#include "statmap.h"
#include "mapped_read.h"
#include "find_candidate_mappings.h"
#include "index_genome.h"
#include "quality.h"
#include "iterative_mapping.h"
#include "candidate_mapping.h"
#include "fragment_length.h"
#include "sam.h"
#include "wiggle.h"
#include "trace.h"
#include "pseudo_location.h"
#include "diploid_map_data.h"
#include "util.h"
#include "error_correction.h"
#include "r_lib.h"

/* Set "unset" defaults for these two global variables */
int num_threads = -1;
int min_num_hq_bps = -1;

/* Getters and setters for utilities written using ctypes */
int get_num_threads() { return num_threads;}
void set_num_threads(int n) { num_threads = n; }

/*
 * Loads the binary genome file
 */
void load_genome( struct genome_data** genome, struct args_t* args )
{
    int rv;

    /* test for a correctly converted binary genome */
    FILE* genome_fp = fopen( args->genome_fname, "r" );
    if( NULL == genome_fp )
    {
        fprintf(stderr, "FATAL       :  Unable to open genome file '%s'\n", args->genome_fname);
        exit( 1 );
    }

    char magic_number[9];
    rv = fread( magic_number, sizeof(char), 9, genome_fp );
    assert( rv == 9 );
    fclose( genome_fp );

    if( 0 == strncmp( magic_number, "SM_OD_GEN", 9 ) )
    {
        /* copy the genome into the output directory */
        safe_link_into_output_directory( 
            args->genome_fname, "./", GENOME_FNAME );
        /* copy genome index into the output directory */
        safe_link_into_output_directory( 
            args->genome_index_fname, "./", GENOME_INDEX_FNAME );
        
        /* copy pseudo_locs into the output directory */
        char pseudo_loc_ofname[500];
        sprintf(pseudo_loc_ofname, "%s.pslocs", args->genome_index_fname  );
        safe_link_into_output_directory( 
            pseudo_loc_ofname, "./", GENOME_INDEX_PSLOCS_FNAME );

        /* copy diploid map data into the output directory */
        char dmap_ofname[500];
        sprintf( dmap_ofname, "%s.dmap", args->genome_index_fname );
        safe_link_into_output_directory(
            dmap_ofname, "./", GENOME_INDEX_DIPLOID_MAP_FNAME );
        
        load_genome_from_disk( genome, GENOME_FNAME );
    }
    else
    {
        fprintf(stderr, "FATAL      :  Genome file '%s' is not in the correct format.\n", args->genome_fname);
        exit(1);
    }

    load_ondisk_index( GENOME_INDEX_FNAME, &((*genome)->index) );

    return;
}

void
map_marginal( struct args_t* args, 
              struct genome_data* genome, 
              struct rawread_db_t* rdb,
              struct mapped_reads_db** mpd_rds_db,
              enum bool is_nc )
{
    /* store clock times - useful for benchmarking */
    struct timeval start, stop;
    
    /* 
       if the error data is not initalized, then we need to bootstrap it. We 
       do this by mapping the reads using a mismatch procedure until we have 
       enough to estiamte the errors
    */
    // bootstrap_error_data

    struct error_model_t* error_model = NULL;
    if( args->search_type == ESTIMATE_ERROR_MODEL )
    {
        fprintf(stderr, "NOTICE      :  Bootstrapping error model\n" );
        init_error_model( &error_model, ESTIMATED );
        bootstrap_estimated_error_model( 
            genome,
            rdb,
            *mpd_rds_db,
            error_model
        );
        
        /* rewind rawread db to beginning for mapping */
        rewind_rawread_db( rdb );
    } else if(args->search_type == PROVIDED_ERROR_MODEL) {
        fprintf(stderr, "FATAL       :  PROVIDED_ERROR_MODEL is not implemented yet\n" );
        exit( 1 );
    } else if(args->search_type == MISMATCHES) {
        init_error_model( &error_model, MISMATCH );
    } else {
        fprintf(stderr, "FATAL       :  Unrecognized index search type '%i'\n",
            args->search_type);
        exit( 1 );
    }
    
    /* initialize the mapped reads db */
    if( false == is_nc )
    {
        open_mapped_reads_db_for_writing( mpd_rds_db, MAPPED_READS_DB_FNAME );
    } else {
        open_mapped_reads_db_for_writing( mpd_rds_db, MAPPED_NC_READS_DB_FNAME );
    }
    
    fprintf(stderr, "NOTICE      :  Finding candidate mappings.\n" );    

    find_all_candidate_mappings( 
        genome,
        rdb,
        *mpd_rds_db,
        error_model,
        args->min_match_penalty,
        args->max_penalty_spread
    );
    
    free_error_model( error_model );

    /* close the mapped reads db that was opened for writing, and reopen it for
     * reading */
    close_mapped_reads_db( mpd_rds_db );

    if( false == is_nc )
    {
        open_mapped_reads_db_for_reading( mpd_rds_db, MAPPED_READS_DB_FNAME );
    } else {
        open_mapped_reads_db_for_reading( mpd_rds_db, MAPPED_NC_READS_DB_FNAME );
    }

    /* write the non-mapping reads into their own fastq */
    gettimeofday( &start, NULL );
    fprintf(stderr, "NOTICE      :  Writing non mapping reads to FASTQ files.\n" );
    write_nonmapping_reads_to_fastq( rdb, *mpd_rds_db );
    gettimeofday( &stop, NULL );
    fprintf(stderr, "PERFORMANCE :  Wrote non-mapping reads to FASTQ in %.2lf sec\n", 
                    (float)(stop.tv_sec - start.tv_sec) 
                        + ((float)(stop.tv_usec - start.tv_usec))/1000000 );

    return;
}

void
build_fl_dist( struct args_t* args, struct mapped_reads_db* mpd_rds_db )
{
    /* estimate the fragment length distribution */
    /* if necessary. and if there are paired reads */
    if( args->frag_len_fp == NULL 
        && args->pair1_reads_fnames != NULL )
    {
        fprintf(stderr, "NOTICE      :  Estimating FL distribution\n" );
        estimate_fl_dist_from_mapped_reads( mpd_rds_db );
    }

    /* write the estimated FL distribution to file */
    if( NULL != mpd_rds_db->fl_dist )
    {
        FILE* fp = fopen( "estimated_fl_dist.txt", "w" );
        if( fp == NULL )
        {
            perror( "ERROR       :  Can not open 'estimated_fl_dist.txt' for writing " );
        } else {
            fprint_fl_dist( fp, mpd_rds_db->fl_dist );
            fclose( fp );
        }
    }
    
    return;
}

void
iterative_mapping( struct args_t* args, 
                   struct genome_data* genome, 
                   struct mapped_reads_db* mpd_rds_db )
{   
    if( NULL != args->unpaired_reads_fnames 
        && args->assay_type == CHIP_SEQ
        && mpd_rds_db->fl_dist == NULL )
    {
        fprintf(stderr, "FATAL       :  Can not iteratively map single end chip-seq reads unless a FL dist is provided\n");
        exit(-1);
    }
 
    /* Do the iterative mapping */
    generic_update_mapping( mpd_rds_db,
                            genome, 
                            args->assay_type,
                            args->num_starting_locations, 
                            MAX_PRB_CHANGE_FOR_CONVERGENCE );
        
    return;
}

void
map_generic_data(  struct args_t* args )
{
    struct genome_data* genome;
    
    /* store clock times - useful for benchmarking */
    struct timeval start, stop;    
    
    /***** Index the genome */
    gettimeofday( &start, NULL );
    
    /* initialize the genome */
    load_genome( &genome, args );
    
    gettimeofday( &stop, NULL );
    fprintf(stderr, "PERFORMANCE :  Indexed Genome in %.2lf seconds\n", 
                    (float)(stop.tv_sec - start.tv_sec) 
                        + ((float)(stop.tv_usec - start.tv_usec))/1000000 );
        
    /***** END Genome processing */
    
    struct mapped_reads_db* mpd_rds_db;    
    map_marginal( args, genome, args->rdb, &mpd_rds_db, false );
    
    /* Free the genome index */
    /* we may need the memory later */
    fprintf(stderr, "NOTICE      :  Freeing index\n" );
    free_ondisk_index( genome->index );
    genome->index = NULL;
    
    /* iterative mapping */
    iterative_mapping( args, genome, mpd_rds_db );
    
    if( args->frag_len_fp != NULL ) {
        build_fl_dist_from_file( mpd_rds_db, args->frag_len_fp );
    } else {
        build_fl_dist( args, mpd_rds_db );
    }

    close_mapped_reads_db( &mpd_rds_db );
    
    free_genome( genome );
    
    return;
}

void
map_chipseq_data(  struct args_t* args )
{
    struct genome_data* genome;

    /* store clock times - useful for benchmarking */
    struct timeval start, stop;    
    
    /***** Index the genome */
    gettimeofday( &start, NULL );
    /* initialize the genome */
    load_genome( &genome, args );
    gettimeofday( &stop, NULL );
    
    fprintf(stderr, "PERFORMANCE :  Indexed Genome in %.2lf seconds\n", 
                    (float)(stop.tv_sec - start.tv_sec) 
                        + ((float)(stop.tv_usec - start.tv_usec))/1000000 );
        
    /***** END Genome processing */
    
    /* map the real ( IP ) data */
    struct mapped_reads_db* chip_mpd_rds_db = NULL;    
    map_marginal( args, genome, args->rdb, &chip_mpd_rds_db, false );
    
    if( args->frag_len_fp != NULL ) {
        build_fl_dist_from_file( chip_mpd_rds_db, args->frag_len_fp );
    } else {
        build_fl_dist( args, chip_mpd_rds_db );
    }
    
    /* 
       this is a bit hacky - for single end chipseq we need to 
       do a bit of work in advance to speed up the fragment 
       coverage smoothing. We do that in the next line.
    */
    if( args->unpaired_reads_fnames != NULL )
        build_chipseq_bs_density( chip_mpd_rds_db->fl_dist );

    struct mapped_reads_db* NC_mpd_rds_db = NULL;
    if ( args->NC_rdb != NULL )
    {        
        map_marginal( args, genome, args->NC_rdb, &NC_mpd_rds_db, true );
        
        if( args->frag_len_fp != NULL ) {
            build_fl_dist_from_file( NC_mpd_rds_db, args->frag_len_fp );
        } else {
            build_fl_dist( args, NC_mpd_rds_db );
        }

        /* 
           this is a bit messy - for single end chipseq we need to 
           do a bit of work in advance to speed up the fragment 
           coverage smoothing. We do that in the next line.
        */
        if( args->unpaired_reads_fnames != NULL )
            build_chipseq_bs_density( NC_mpd_rds_db->fl_dist );
    }

    /* Free the genome index */
    /* we may need the memory later */
    fprintf(stderr, "NOTICE      :  Freeing index\n" );
    free_ondisk_index( genome->index );
    
    /* if there is no negative control, we use the same iterative 
       scheme as the generic version. iterative_mapping takes care of 
       everything ( output, iterative, etc. ). We dont touch peak calling - 
       ( we dont really know how to do it well inside our probability model
         without a NC because of chromatin solubility, etc., so we leave that for
         people that have taken the time to build effective heiristics. ie. 
         peak callers.  )
    */
    if( NULL == args->NC_rdb )
    {
        iterative_mapping( args, genome, chip_mpd_rds_db );
    } 
    /* if there is a NC, we can do something simple. For each 
       sample that we take from the mapping distribution, we use the 
       machinery to build a marginal read density. Then, we update the 
       read mapping locations for the NC control from the marginal density, 
       and build a NC marginal density. Of course, the marginal density is 
       precisely the expectation of the binding site density distribution,
       so calling peaks ( basepair per basepair ) corresponds to testing the
       hypothesis that the binding site in the IP sample is greater than 
       the negative control. 
    */
    else {
        /* Iterative mapping */
        /* mmap and index the necessary data */
        fprintf(stderr, "NOTICE      :  mmapping mapped IP reads DB.\n" );
        mmap_mapped_reads_db( chip_mpd_rds_db );
        fprintf(stderr, "NOTICE      :  indexing mapped IP reads DB.\n" );
        index_mapped_reads_db( chip_mpd_rds_db );

        fprintf(stderr, "NOTICE      :  mmapping mapped NC reads DB.\n" );
        mmap_mapped_reads_db( NC_mpd_rds_db );
        fprintf(stderr, "NOTICE      :  indexing mapped NC reads DB.\n" );
        index_mapped_reads_db( NC_mpd_rds_db );
        
        /* open the file to store the meta info, and print out the header */
        FILE* s_mi = fopen(RELAXED_SAMPLES_META_INFO_FNAME, "w");
        fprintf( s_mi, "sample_number,log_lhd\n" );
        
        int i;
        for( i = 0; i < args->num_starting_locations; i++ )
        {
            take_chipseq_sample_wnc(
                chip_mpd_rds_db, NC_mpd_rds_db,
                genome, 
                s_mi,  // meta info file pointer
                i, // sample index
                MAX_PRB_CHANGE_FOR_CONVERGENCE, 
                true // use a random starting location
            );
        }
    
        fclose( s_mi );
    }
    
    goto cleanup;

cleanup:    
    close_mapped_reads_db( &chip_mpd_rds_db );
    close_mapped_reads_db( &NC_mpd_rds_db );

    free_genome( genome );
    
    return;
}

int 
main( int argc, char** argv )
{       
    /* Seed the random number generator */
    srand ( time(NULL) );

    /* get the base directory */
    char* abs_path = NULL;
    abs_path = realpath( argv[0], abs_path );
    if( NULL == abs_path )
    {
        fprintf( stderr, "%s\n", argv[0] );
        perror( "Couldnt find abs path" );
        exit( -1 );
    }
    char* statmap_base_dir = dirname( abs_path );

    /* parse and sanity check arguments */
    struct args_t args = parse_arguments( argc, argv );
    
    /* write the arguments to a file */
    FILE* arg_fp = fopen( "config.dat", "w" );
    write_config_file_to_stream( arg_fp, &args );
    fclose( arg_fp  );
    
    /* intialize an R instance */
    if( args.search_type == ESTIMATE_ERROR_MODEL )
    {
        fprintf( stderr, "NOTICE      :  Initializing R interpreter.\n" );
        init_R( );        
        load_statmap_source( statmap_base_dir );
    }
    
    if( args.assay_type == CHIP_SEQ )
    {
        map_chipseq_data( &args );
    } else {
        map_generic_data( &args );
    }
       
    goto cleanup;
    
cleanup:
    close_rawread_db( args.rdb );
    
    if( args.NC_rdb != NULL )
        close_rawread_db( args.NC_rdb );
    
    if( args.frag_len_fp != NULL ) {
        fclose( args.frag_len_fp );
    }

    free( args.genome_fname );
    free( args.genome_index_fname );
    free( args.output_directory );
    
    /* finish the R interpreter */
    end_R();
    fprintf( stderr, "NOTICE      :  Shutting down R interpreter.\n" );

    return 0;
}
