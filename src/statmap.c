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

#include "igraph.h"

#include "log.h"
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
int max_reference_insert_len = -1;

/* TODO check for conflicts with variable name before removing 
   leading underscore */
enum assay_type_t _assay_type = UNKNOWN;
/* TODO should this be global? */
int softclip_len = 0;

/* Getters and setters for utilities written using ctypes */
int get_num_threads() { return num_threads;}
void set_num_threads(int n) { num_threads = n; }
int get_max_reference_insert_len() { return max_reference_insert_len; }
void set_max_reference_insert_len(int n) { max_reference_insert_len = n; }

/*
 * Loads the binary genome file
 */
void load_genome( struct genome_data** genome, struct args_t* args )
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

    /* Save the user-specified mapping metaparameters in a struct */
    struct mapping_metaparams mapping_metaparams;

    /* 
       if the error data is not initalized, then we need to bootstrap it. We 
       do this by mapping the reads using a mismatch procedure until we have 
       enough to estiamte the errors
    */
    struct error_model_t* error_model = NULL;
    if( args->error_model_type == ESTIMATED )
    {
        statmap_log( LOG_NOTICE, "Bootstrapping error model" );
        init_error_model( &error_model, ESTIMATED );
        
        statmap_log( LOG_NOTICE,"Setting bootstrap mismatch rates to %f and %f",
                MAX_NUM_MM_RATE, MAX_NUM_MM_SPREAD_RATE  );

        mapping_metaparams.error_model_type = MISMATCH;
        mapping_metaparams.error_model_params[0] = MAX_NUM_MM_RATE;
        mapping_metaparams.error_model_params[1] = MAX_NUM_MM_SPREAD_RATE;
        
        bootstrap_estimated_error_model( 
            genome,
            rdb,
            &mapping_metaparams,
            error_model
        );

        mapping_metaparams.error_model_type = ESTIMATED;
        mapping_metaparams.error_model_params[0] = args->mapping_metaparameter;
        
        /* rewind rawread db to beginning for mapping */
        rewind_rawread_db( rdb );
    } else if(args->error_model_type == FASTQ_MODEL) {
        statmap_log( LOG_FATAL, 
            "FASTQ_MODEL (provided error scores) is not implemented yet" );
        exit( 1 );
    } else if(args->error_model_type == MISMATCH) {
        /* initialize the meta params */
        mapping_metaparams.error_model_type = MISMATCH;
        mapping_metaparams.error_model_params[0]=args->mapping_metaparameter;
        mapping_metaparams.error_model_params[1]=args->mapping_metaparameter/2;
        
        init_error_model( &error_model, MISMATCH );
    } else {
        statmap_log( LOG_FATAL, "Invalid index search type '%i'",  
                     args->error_model_type );
        assert( false );
        exit( 1 );
    }

    /* initialize the mapped reads db */
    if( false == is_nc ) {
        open_mapped_reads_db_for_writing( mpd_rds_db, MAPPED_READS_DB_FNAME );
    } else {
        open_mapped_reads_db_for_writing( mpd_rds_db, MAPPED_NC_READS_DB_FNAME );
    }

    /* initialize the fl dist if one was provided. This is necessary to properly
       compute the penalty for a read, including the penalty from the fragment
       length */
    if( args->frag_len_fp != NULL ) {
        init_fl_dist_from_file( &((*mpd_rds_db)->fl_dist), args->frag_len_fp );
    }
    
    statmap_log( LOG_NOTICE, "Finding candidate mappings." );
    find_all_candidate_mappings( 
        genome,
        rdb,
        *mpd_rds_db,
        &mapping_metaparams,
        error_model
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

    /* Reload the fragment length distribution
       TODO: there might be a cleaner way to do all of this */
    if( args->frag_len_fp != NULL ) {
        init_fl_dist_from_file( &((*mpd_rds_db)->fl_dist), args->frag_len_fp );
    }

    /* write the non-mapping reads into their own fastq */
    gettimeofday( &start, NULL );
    statmap_log( LOG_NOTICE, "Writing non mapping reads to FASTQ files." );
    write_nonmapping_reads_to_fastq( rdb, *mpd_rds_db );
    gettimeofday( &stop, NULL );
    statmap_log( LOG_INFO, "Wrote non-mapping reads to FASTQ in %.2lf sec",
            (float)(stop.tv_sec - start.tv_sec) + ((float)(stop.tv_usec - start.tv_usec))/1000000 );

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
        statmap_log( LOG_NOTICE, "Estimating FL distribution" );
        estimate_fl_dist_from_mapped_reads( mpd_rds_db );
    }

    /* write the estimated FL distribution to file */
    if( NULL != mpd_rds_db->fl_dist )
    {
        FILE* fp = fopen( "estimated_fl_dist.txt", "w" );
        if( fp == NULL )
        {
            statmap_log( LOG_ERROR, "Can not open 'estimated_fl_dist.txt' for writing" );
        } else {
            fprint_fl_dist( fp, mpd_rds_db->fl_dist );
            fclose( fp );
        }
    }
    
    return;
}

int 
main( int argc, char** argv )
{
    /* Seed the random number generator */
    srand ( time(NULL) );

    /* store clock times - useful for benchmarking */
    struct timeval start, stop;
    
    /* Turn on attribute handling.
       Manual recommends "set at the beginning of main and never touch again":
       http://igraph.sourceforge.net/doc/html/ch12s02.html */
    igraph_i_set_attribute_table(&igraph_cattribute_table);

    /* Log to a temporary log file until we are ready to create the real log
     * file (must wait until the Statmap output directory has been created) */
    init_initial_logging();

    /* get the base directory */
    char* abs_path = NULL;
    abs_path = realpath( argv[0], abs_path );
    if( NULL == abs_path )
    {
        statmap_log( LOG_FATAL, "Couldn't find abs path: %s", argv[0] );
    }
    char* statmap_base_dir = dirname( abs_path );

    /* parse and sanity check arguments */
    struct args_t args = parse_arguments( argc, argv );
    
    /* write the arguments to a file */
    FILE* arg_fp = fopen( "config.dat", "w" );
    write_config_file_to_stream( arg_fp, &args );
    fclose( arg_fp  );
    
    /* intialize an R instance */
    if( args.error_model_type == ESTIMATED )
    {
        statmap_log( LOG_NOTICE, "Initializing R interpreter." );
        init_R( );        
        load_statmap_source( statmap_base_dir );
    }

    struct genome_data* genome;
        
    /* load the genome */
    gettimeofday( &start, NULL );
    load_genome( &genome, &args );
    gettimeofday( &stop, NULL );
    statmap_log( LOG_INFO, "Loaded Genome in %.2lf seconds",
            (float)(stop.tv_sec - start.tv_sec) 
                 + ((float)(stop.tv_usec - start.tv_usec))/1000000 );
    
    /***** END Genome processing */
    
    struct mapped_reads_db* mpd_rds_db;    
    map_marginal( &args, genome, args.rdb, &mpd_rds_db, false );

    /* Free the genome index */
    /* we may need the memory later */
    free_ondisk_index( genome->index );
    genome->index = NULL;
    
    /* iterative mapping */
    if( args.error_model_type == ESTIMATED
        && args.assay_type != UNKNOWN )
    {
        statmap_log( LOG_NOTICE, "Starting iterative mapping" );
        struct cond_prbs_db_t* cond_prbs_db = build_posterior_db( 
            genome, mpd_rds_db, args.assay_type );
    
        // args.sam_output_fname
        statmap_log( LOG_NOTICE, "Writing mapped reads to sam file" );
        FILE* fp = fopen("mapped_reads.sam", "w");
        write_mapped_reads_to_sam( 
            args.rdb, mpd_rds_db, cond_prbs_db, genome, false, fp );
        fclose(fp);

        free_cond_prbs_db(cond_prbs_db);
    } else {
        statmap_log( LOG_NOTICE, "Writing mapped reads to sam file" );
        struct cond_prbs_db_t* cond_prbs_db;
        init_cond_prbs_db_from_mpd_rdb( &cond_prbs_db, mpd_rds_db);
        FILE* fp = fopen("mapped_reads.sam", "w");
        write_mapped_reads_to_sam( 
            args.rdb, mpd_rds_db, cond_prbs_db, genome, true, fp );
        free_cond_prbs_db(cond_prbs_db);
        fclose(fp);

    }
    goto cleanup;
    
cleanup:
    free( abs_path );

    close_mapped_reads_db( &mpd_rds_db );
    close_rawread_db( args.rdb );    
    free_genome( genome );

    if( args.NC_rdb != NULL )
        close_rawread_db( args.NC_rdb );
    
    if( args.frag_len_fp != NULL ) {
        fclose( args.frag_len_fp );
    }

    free( args.genome_fname );
    free( args.genome_index_fname );
    free( args.output_directory );
    
    /* finish the R interpreter */
    statmap_log( LOG_NOTICE, "Shutting down R interpreter." );
    end_R();

    finish_logging();

    return 0;
}
