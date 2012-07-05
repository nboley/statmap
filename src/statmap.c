/* Copyright (c) 2009-2010, Nathan Boley */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h> /* gettimeofday() */

#include <getopt.h>
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
#include "error_correction.h"

/* Set the deafults for these two global variables */
int num_threads = -1;
int min_num_hq_bps = -1;

/* Getters and setters for utilities written using ctypes */
int get_num_threads() { return num_threads;}
void set_num_threads(int n) { num_threads = n; }
int get_min_num_hq_bps() { return min_num_hq_bps; }
void set_min_num_hq_bps(int n) { min_num_hq_bps = n; }

// TODO: update message
void usage()
{
    fprintf(stderr, "Usage: ./statmap -g genome.fa.bin -p men_match_penalty -m max_penalty_spread \n");
    fprintf(stderr, "                 ( -r input.fastq | [ -1 input.pair1.fastq & -2 input.pair2.fastq ] ) \n");
    fprintf(stderr, "      (optional)  [ -o output_directory  -a assay_type -f fragment_lengths ] \n\n" );
}

void
safe_mkdir(char* dir)
{
    int error = mkdir( dir, S_IRWXU | S_IRWXG | S_IRWXO );
    char buffer[200];
    sprintf( buffer, "FATAL       :  Cannot make %s ", dir );
    if( -1 == error )
    {
        perror( buffer );
        assert( false );
        exit( -1 );
    }
}

static void
safe_copy_into_output_directory( char* fname, char* output_dir, char* output_fname )
{
    int error;

    char* realpath_buffer = realpath( output_dir, NULL );
    assert( NULL != realpath_buffer );

    char buffer[500];
    sprintf( buffer, "cp %s %s/%s", fname, realpath_buffer, output_fname );
    fprintf(stderr, "NOTICE      :  Copying '%s' to the output directory\n",  fname );
    error = system( buffer );
    if (WIFSIGNALED(error) &&
        (WTERMSIG(error) == SIGINT || WTERMSIG(error) == SIGQUIT))
    {
        fprintf(stderr, "FATAL     : Failed to call '%s'\n", buffer );
        perror( "Copy failure");
        assert( false );
        exit( -1 );
    }
    
    free( realpath_buffer );
    
    return;
}

static void
safe_link_into_output_directory( char* fname, char* output_dir, char* output_fname )
{
    int error;

    char* realpath_buffer = realpath( output_dir, NULL );
    assert( NULL != realpath_buffer );

    char* input_realpath_buffer = realpath( fname, NULL );
    assert( NULL != realpath_buffer );

    char buffer[500];
    sprintf( buffer, "ln -s %s %s/%s", input_realpath_buffer, realpath_buffer, output_fname );
    fprintf(stderr, "NOTICE      :  Linking '%s' to the output directory\n",  fname );
    error = system( buffer );
    if (WIFSIGNALED(error) &&
        (WTERMSIG(error) == SIGINT || WTERMSIG(error) == SIGQUIT))
    {
        fprintf(stderr, "FATAL     : Failed to call '%s'\n", buffer );
        perror( "Link failure");
        assert( false );
        exit( -1 );
    }
    
    free( realpath_buffer );
    free( input_realpath_buffer );
    
    return;
}

enum input_file_type_t
guess_input_file_type( struct args_t* args )
{
    enum input_file_type_t input_file_type = UNKNOWN;

    // Store the minimum and max qual scores among
    // every observed basepair.
    unsigned char max_qual = 0;
    unsigned char min_qual = 255;

    // Open one of of the read input files
    char* if_name = args->unpaired_reads_fnames;
    if( if_name == NULL )
        if_name = args->pair1_reads_fnames;
    FILE* fp = fopen( if_name, "r" );

    /* FIXME - stop assuming this is a fastq file */
    int i, j;
    for( i = 0; i < 10000; i++ )
    {
        struct rawread* r;
        populate_read_from_fastq_file( fp, &r, UNKNOWN );
        /* If we have reached and EOF, then break */
        if( r == NULL )
            break;

        for( j = 0; j < r->length; j++ )
        {
            max_qual = MAX( max_qual, (unsigned char ) r->error_str[j]  );
            min_qual = MIN( min_qual, (unsigned char ) r->error_str[j]  );
        }
        free_rawread( r );
    }

    fprintf( stderr,
             "NOTICE      :  Calculated max and min quality scores as '%i' and '%i'\n",
             max_qual, min_qual );


    /* if the max quality score is 73, ( 'I' ), then
       this is almost certainly a sanger format */
    if( max_qual == 73 ) {
        input_file_type = SANGER_FQ ;
    } else if( max_qual > 73 && max_qual <= 90 )
    {
        input_file_type = MARKS_SOLEXA;
    } else {
        // If the max quality is 104 it's a bit tougher.
        // because it can be either the new or old illumina format
        if( max_qual > 73 ) {
            if( max_qual < 104 )
            {
                /* If the maximum quality is less than 104 but greater than 73, it's probably the
                   plus 64 version, but we can't be sure. However, assume that it is and print a 
                   warning. 
                */
                fprintf(stderr, "WARNING     :  maximum input score wasn't achieved. ( %i )\n", max_qual);
                fprintf(stderr, "WARNING     :  ( That means the highest quality basepair was above 73 but less than 104. Probably, there are just no very HQ basepairs, but make sure that the predicted error format is correct. )\n");
            }
            /*
             * If the shifted min quality is less than 0, then the
             * format is almost certainly the log odds solexa version
             */
            if( min_qual < 64 )
            {
                input_file_type = SOLEXA_LOG_ODDS_FQ;
            }
            // If the min is less than 69, it is still most likely
            // a log odds version, but we emit a warning anyways.
            else if ( min_qual < 66 ) {
                fprintf(stderr, "WARNING     :  input file format ambiguous. ( max=%i, min=%i )\n", max_qual, min_qual);
                input_file_type = ILLUMINA_v13_FQ;
            } else if ( min_qual < 70 ) {
                fprintf(stderr, "WARNING     :  input file format ambiguous. ( max=%i, min=%i )\n", max_qual, min_qual);
                input_file_type = ILLUMINA_v15_FQ;
            } else if ( min_qual < 90 ) {
                fprintf(stderr, "WARNING     :  input file format ambiguous. ( max=%i, min=%i )\n", max_qual, min_qual);
                input_file_type = MARKS_SOLEXA;
            } else if ( min_qual == 104 ) {
                fprintf(stderr, "WARNING     :  input file format ambiguous. ( max=%i, min=%i )\n", max_qual, min_qual);
                input_file_type = TEST_SUITE_FORMAT;
            } else {
                fprintf( stderr, "ERROR       :  Could not automatically determine " );
                fprintf( stderr, "input format. ( max=%i, min=%i  )\n",
                         max_qual, min_qual );
                exit( -1 );
            }
        }
    }

    /* print out the determined format */
    fprintf( stderr,
             "NOTICE      :  Setting Input File Format to '%i'\n",
             input_file_type );
    fclose( fp );

    switch( input_file_type )
    {
    case 1:
        QUAL_SHIFT = 33;
        ARE_LOG_ODDS = false;
        break;
    case 2:
        QUAL_SHIFT = 64;
        ARE_LOG_ODDS = false;
        break;
    case 3:
        QUAL_SHIFT = 59;
        ARE_LOG_ODDS = true;
        break;
    case 4:
        QUAL_SHIFT = 50;
        ARE_LOG_ODDS = false;
        break;
    default:
        fprintf( stderr,
                 "ERROR       :  Unrecognized file format type %i\n",
                 input_file_type );
        exit( -1 );
    }

    fprintf( stderr,
             "NOTICE      :  Set Qual_shift to '%i' and ARE_LOG_ODDS to %i\n",
             QUAL_SHIFT, ARE_LOG_ODDS );


    return input_file_type;
}


int
guess_optimal_indexed_seq_len( struct args_t* args)
{
    int seq_len = -1;

    char* if_name = args->unpaired_reads_fnames;
    if( if_name == NULL )
        if_name = args->pair1_reads_fnames;
    
    FILE* fp = fopen( if_name, "r" );
    if( fp == NULL ) {
        fprintf(stderr, "FATAL       :  Unable to open reads file '%s'\n", if_name);
        exit(-1);
    }
    
    /* FIXME - stop assuming this is a fastq file */
    /* TODO - read in multiple reads to corroborate read length */
    /* read in the first read, and set the seq length from this read */
    struct rawread* r;
    populate_read_from_fastq_file( fp, &r, UNKNOWN );
    seq_len = r->length;
    fprintf( stderr, 
             "NOTICE      :  Setting Read Length to '%i' BPS\n", 
             seq_len );
    
    free_rawread( r );
    fclose( fp );
    
    return seq_len;
}

struct args_t
parse_arguments( int argc, char** argv )
{
    int error;

    /* 
     * initialize the structure that stores all of the 
     * configuration options.
     */
    struct args_t args;
    
    args.genome_fname = NULL;
    args.genome_index_fname = NULL;
    
    args.unpaired_reads_fnames = NULL;
    args.pair1_reads_fnames = NULL;
    args.pair2_reads_fnames = NULL;
    args.rdb = NULL;

    args.unpaired_NC_reads_fnames = NULL;
    args.pair1_NC_reads_fnames = NULL;
    args.pair2_NC_reads_fnames = NULL;
    args.NC_rdb = NULL;

    args.frag_len_fname = NULL;
    args.frag_len_fp = NULL;
    
    args.output_directory = NULL;
    
    args.sam_output_fname = NULL;
    
    args.min_match_penalty = 1;
    args.max_penalty_spread = -1;
    args.min_num_hq_bps = -1;

    args.num_starting_locations = DEFAULT_NUM_SAMPLES;

    args.num_threads = -1;
    
    args.input_file_type = UNKNOWN;
    args.assay_type = UNKNOWN;

    char* assay_name = NULL;
    char* search_type = NULL;

    int c;
    while ((c = getopt(argc, argv, "9hg:n:r:1:2:3:4:c:o:sp:m:f:l:a:w:t:q:y:")) != -1) {
        switch (c) {
        /* Required Argumnets */
        case 'g': // reference genome fasta file
            args.genome_fname = optarg;
            break;
        case 'p': // minimum match penalty
            args.min_match_penalty = atof(optarg);
            break;
        case 'm': // maximum penalty spread
            args.max_penalty_spread = atof(optarg);
            break;
        
        /* input file arguments */
        /* either ( -1 and -2 ) or -r is required */
        case '1': // paired end reads input file 1
            args.pair1_reads_fnames = optarg;
            break;
        case '2': // paired end reads input file 2
            args.pair2_reads_fnames = optarg;
            break;
        case 'r': // single end reads input file
            args.unpaired_reads_fnames = optarg;
            break;

        /* either ( -3 and -4 ) or -c should be set exclusively, but neither are required */
        case '3': // paired end reads input file 1
            args.pair1_NC_reads_fnames = optarg;
            break;
        case '4': // paired end reads input file 2
            args.pair2_NC_reads_fnames = optarg;
            break;
        case 'c': // single end reads input file
            args.unpaired_NC_reads_fnames = optarg;
            break;
            
        /* optional arguments ( that you probably want ) */
        case 'o': // output directory
            // since -o might be dynamically allocated (if this argument is not set),
            // we copy the string here so we can properly free memory when we're done
            args.output_directory = calloc( strlen(optarg)+1, sizeof( char ) );
            strcpy( args.output_directory, optarg );
            break;
        case 's': // whether to write sam files or not
            args.sam_output_fname = SAM_MARGINAL_OFNAME;
            break;
        
        case 'a': // the assay type
            assay_name = optarg;
            break;

        case 'f': // fragment length file
            args.frag_len_fname = optarg;
            break;
        
        case 'q': // min number of HQ basepairs
            args.min_num_hq_bps = atoi( optarg );
            break;
        case 'n': // number of starting locations
            args.num_starting_locations = atoi(optarg);
            break;
        
        case 't': // number of threads
            args.num_threads = atoi( optarg );
            break;

        case 'y': // type of marginal mapping to do - mismatches, or error data
            search_type = optarg;
            break;
        
        /* utility options */
        case 'h':
            usage();
            exit(-1);            
        case '?':
            fprintf(stderr, "FATAL       :  Unrecognized Argument: '%c' \n", (char) optopt);
            usage();
            exit(-1);
        default:
            usage();
            exit( -1 );
        }
    }

    /********* CHECK REQUIRED ARGUMENTS *************************************/
    /* Ensure that the required arguments were present */

    if( args.genome_fname == NULL )
    {
        usage();
        fprintf(stderr, "FATAL       :  -g ( binary genome ) is required\n");
        exit( -1 );
    }
    
    if( args.unpaired_reads_fnames == NULL
        && ( args.pair2_reads_fnames == NULL
             || args.pair1_reads_fnames == NULL ) 
        )
    {
        usage();
        fprintf(stderr, "FATAL       :  -r or ( -1 and -2 ) is required\n" );
        exit( -1 );
    }
    
    /* perform sanity checks on the read input arguments */
    if( args.unpaired_reads_fnames )
    {
        if( args.pair1_reads_fnames != NULL
            || args.pair2_reads_fnames != NULL     )
        {
            usage();
            fprintf(stderr, 
             "FATAL       :  if the -r is set, neither -1 nor -2 should be set\n"
            );
            exit(-1);
        }
    }

    if( args.pair1_reads_fnames != NULL 
        && args.pair2_reads_fnames == NULL )
    {
        usage();
        fprintf(stderr, 
                "FATAL       :  -1 option is set but -2 is not\n" );
        exit(-1);            
    }

    if( args.pair2_reads_fnames != NULL 
        && args.pair1_reads_fnames == NULL )
    {
        usage();
        fprintf(stderr, 
                "FATAL       :  -2 option is set but -1 is not\n" );
        exit(-1);            
    }

    /********* END CHECK REQUIRED ARGUMENTS ************************************/

    /* set the assay type  */
    if( assay_name != NULL )
    {
        switch (assay_name[0] )
        {
        case 'a':
            args.assay_type = CAGE;
            break;
        case 'i':
            args.assay_type = CHIP_SEQ;
            break;
        default:
            fprintf(stderr, 
                    "FATAL       :  Unrecognized Assay Type: '%s'\n", 
                    assay_name );
            assert( false );
            exit( -1 );
        }
    }

    /* set the search type */
    if( search_type != NULL )
    {
        switch( search_type[0] )
        {
        case 'e':
            args.search_type = ESTIMATE_ERROR_MODEL;
            break;
        case 'm':
            args.search_type = MISMATCHES;
            break;
        default:
            fprintf(stderr,
                    "FATAL       :  Unrecognized search type: '%s'\n",
                    search_type );
            exit(-1);
        }
    } else {
        /* set default to be using the observed error rates */
        args.search_type = ESTIMATE_ERROR_MODEL;
    }

    /* Make the output directory */
    if( args.output_directory == NULL )
    {
        time_t rawtime;
        time ( &rawtime );
        struct tm * timeinfo;
        timeinfo = localtime ( &rawtime );        
        char buffer[200];
        strftime(buffer, 200, "statmap_output_%Y_%m_%d_%H_%M_%S", timeinfo);
        args.output_directory = calloc(strlen(buffer)+1, sizeof(char));
        memcpy( args.output_directory, buffer, (strlen(buffer)+1)*sizeof(char) );
    } 
    
    error = mkdir( args.output_directory, S_IRWXU | S_IRWXG | S_IRWXO );
    if( -1 == error )
    {
        perror( "FATAL       :  Cannot make output directory");
        exit( -1 );
    }

    /* Check path of genome file
     * Expand paths to absolute paths
     * Set args.genome_index_fname based on args.genome_fname.
     */ 
    args.genome_index_fname = calloc( sizeof(char), PATH_MAX + 6 );
    if( args.genome_fname != NULL ) {
        char* genome_fname = realpath( args.genome_fname, NULL );
        assert( NULL != genome_fname );
        args.genome_fname = genome_fname;
        sprintf( args.genome_index_fname, "%s.index", args.genome_fname );
    }
    
    if( args.min_match_penalty == 1 )
    {
        args.min_match_penalty = DEFAULT_MIN_MATCH_PENALTY;
        fprintf(stderr, "NOTICE      :  Setting the min_match_penalty (-p) to %.3f\n", DEFAULT_MIN_MATCH_PENALTY );        
    }
    
    if( args.max_penalty_spread == -1 )
    {
        args.max_penalty_spread = DEFAULT_MAX_PENALTY_SPREAD;
        fprintf(stderr, "NOTICE      :  Setting the max_penalty_spread (-m) to %.3f\n", DEFAULT_MAX_PENALTY_SPREAD );        
    }
    
    /***** Copy the reads file into the output directory ****/
    /* If the reads are not paired */
    if( args.unpaired_reads_fnames != NULL )
    {
        /* first, copy the read file(s) into the output directory */
        safe_copy_into_output_directory( 
            args.unpaired_reads_fnames, args.output_directory, "reads.unpaired" );

        /* copy the unpaired NC reads */
        if( args.unpaired_NC_reads_fnames != NULL )
        {
            safe_copy_into_output_directory( 
                args.unpaired_NC_reads_fnames, args.output_directory, "reads.NC.unpaired" );
        }
    } 
    /* If the reads are paired */
    else {
        safe_copy_into_output_directory( 
            args.pair1_reads_fnames, args.output_directory, "reads.pair1" );
        
        safe_copy_into_output_directory( 
            args.pair2_reads_fnames, args.output_directory, "reads.pair2" );

        /* if they exist, copy the negatice control reads */ 
        if( args.pair1_NC_reads_fnames != NULL )
        {
            safe_copy_into_output_directory( 
                args.pair1_NC_reads_fnames, args.output_directory, "reads.NC.pair1" );
            
            if( args.pair2_NC_reads_fnames == NULL )
            {
                fprintf( stderr, "FATAL     : The NC reads must be paired for a paired experiment.\n" );
                assert( false );
                exit( 1 );
            }
            
            safe_copy_into_output_directory( 
                args.pair2_NC_reads_fnames, args.output_directory, "reads.NC.pair2" );            
        }
    }

    /*
     * Try to determine the type of sequencing error. 
     */
    if( args.input_file_type == UNKNOWN )
    {
        args.input_file_type = guess_input_file_type( &args );
    }

    
    if( args.frag_len_fname != NULL )
    {
        safe_copy_into_output_directory( 
            args.frag_len_fname, args.output_directory, "estimated_fl_dist.txt" );
        
        args.frag_len_fp = fopen( args.frag_len_fname, "r" );
        if( NULL == args.frag_len_fp )
        {
            char buffer[500];
            sprintf( buffer, "FATAL       :  Failed to open '%s' ", args.frag_len_fname );
            perror( buffer );
            exit( -1 );
        }    
    } else {
        if( args.assay_type == CHIP_SEQ )
        {
            char buffer[500];
            sprintf( buffer, "FATAL       :  CHIPSEQ assay mappings requires a fragment length distirbution (-f)\n" );
            perror( buffer );
            exit( 1 );
        }
    }


    /* change the working directory to the output directory */
    error = chdir( args.output_directory );
    if( -1 == error )
    {
        perror( "FATAL       :  Cannot move into output directory ");
        exit( -1 );
    }

    /* Make sub directories to store wiggle samples */
    if( SAVE_STARTING_SAMPLES )
        safe_mkdir( STARTING_SAMPLES_PATH );
    
    if( SAVE_SAMPLES )
        safe_mkdir( RELAXED_SAMPLES_PATH );

    if( CALL_PEAKS )
        safe_mkdir( CALLED_PEAKS_OUTPUT_DIRECTORY );
    
    if( SAVE_BOOTSTRAP_SAMPLES )
    {
        safe_mkdir( BOOTSTRAP_SAMPLES_PATH );
        safe_mkdir( BOOTSTRAP_SAMPLES_MAX_PATH );
        safe_mkdir( BOOTSTRAP_SAMPLES_MIN_PATH );
        safe_mkdir( BOOTSTRAP_SAMPLES_ALL_PATH );
    }

    /***** initialize the raw reads db */
    init_rawread_db( &(args.rdb) );
    if( args.unpaired_NC_reads_fnames != NULL 
        || args.pair1_NC_reads_fnames != NULL )
    {
        init_rawread_db( &(args.NC_rdb) );
    }
    
    /* If the reads are not paired */
    if( args.unpaired_reads_fnames != NULL )
    {
        add_single_end_reads_to_rawread_db(
            args.rdb, "reads.unpaired", FASTQ, args.assay_type 
        );
        
        if( args.unpaired_NC_reads_fnames != NULL )
        {
            add_single_end_reads_to_rawread_db(
                args.NC_rdb, "reads.NC.unpaired", FASTQ, args.assay_type
            );
        }
    } 
    /* If the reads are paired */
    else {
        add_paired_end_reads_to_rawread_db(
            args.rdb, 
            "reads.pair1",
            "reads.pair2",
            FASTQ,
            args.assay_type
        );

        if( args.pair1_NC_reads_fnames != NULL )
        {
            add_paired_end_reads_to_rawread_db(
                args.NC_rdb, 
                "reads.NC.pair1",
                "reads.NC.pair2",
                FASTQ,
                args.assay_type
            );
        }
    }
    /***** END initialize the read db */


    /* If we didnt set the number of threads, set it to the maximum available */
    if( args.num_threads == -1 )
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
        num_threads = args.num_threads;
    }
    
    /* set the min num hq basepairs if it's unset */
    if( args.min_num_hq_bps == -1 )
    {
        args.min_num_hq_bps = 12;
    }
    min_num_hq_bps = args.min_num_hq_bps;
    
    /* 
     * Dont allow penalty spreads greater than the min match penalty - 
     * they are worthless and mess up my search heuristics 
     *
     */
    if( args.max_penalty_spread + args.min_match_penalty + 0.00001 > 0 )
        args.max_penalty_spread = -1;
    
    return args;
}

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
    
    /***** initialize the mappings dbs */
    
    candidate_mappings_db candidate_mappings;
    if( false == is_nc )
    {
        init_candidate_mappings_db( &candidate_mappings, 
                                    CANDIDATE_MAPPINGS_DB_FNAME );
        new_mapped_reads_db( mpd_rds_db, MAPPED_READS_DB_FNAME );
    } else {
        init_candidate_mappings_db( &candidate_mappings, 
                                    CANDIDATE_MAPPINGS_NC_DB_FNAME );
        new_mapped_reads_db( mpd_rds_db, MAPPED_NC_READS_DB_FNAME );
    }

    
    /***** END initialize the mappings dbs */
    
    /* 
       if the error data is not initalized, then we need to bootstrap it. We 
       do this by mapping the reads using a mismatch procedure until we have 
       enough to estiamte the errors
    */
    // bootstrap_error_data

    struct error_model_t* error_model = NULL;
    if( args->search_type == ESTIMATE_ERROR_MODEL )
    {
        init_error_model( &error_model, ESTIMATED );
        bootstrap_estimated_error_model( 
            genome,
            rdb,
            &candidate_mappings,
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
    
    find_all_candidate_mappings( 
        genome,
        rdb,
        &candidate_mappings,
        error_model,
        args->min_match_penalty,
        args->max_penalty_spread
    );
    
    free_error_model( error_model );
    
    /* combine and output all of the partial mappings - this includes
       joining paired end reads. */
    fprintf(stderr, "NOTICE      :  Joining Candidate Mappings\n" );
    gettimeofday( &start, NULL );
    join_all_candidate_mappings( &candidate_mappings, *mpd_rds_db, genome );
    gettimeofday( &stop, NULL );
    fprintf(stderr, "PERFORMANCE :  Joined Candidate Mappings in %.2lf seconds\n", 
                    (float)(stop.tv_sec - start.tv_sec) 
                        + ((float)(stop.tv_usec - start.tv_usec))/1000000 );
    /*  close candidate mappings db */
    close_candidate_mappings_db( &candidate_mappings );
    
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
 
    /* Iterative mapping */
    /* mmap and index the necessary data */
    fprintf(stderr, "NOTICE      :  mmapping mapped reads DB.\n" );
    mmap_mapped_reads_db( mpd_rds_db );
    fprintf(stderr, "NOTICE      :  indexing mapped reads DB.\n" );
    index_mapped_reads_db( mpd_rds_db );
    
    /* Do the iterative mapping */
    generic_update_mapping( mpd_rds_db,
                            genome, 
                            args->assay_type,
                            args->num_starting_locations, 
                            MAX_PRB_CHANGE_FOR_CONVERGENCE );
        
    munmap_mapped_reads_db( mpd_rds_db );
    
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
    
    /* parse and sanity check arguments */
    struct args_t args = parse_arguments( argc, argv );
    
    /* write the arguments to a file */
    FILE* arg_fp = fopen( "config.dat", "w" );
    write_config_file_to_stream( arg_fp, &args );
    fclose( arg_fp  );
    
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

    return 0;
}
