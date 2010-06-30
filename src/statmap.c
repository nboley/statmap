/* Copyright (c) 2009-2010, Nathan Boley */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <getopt.h>
#include <string.h>
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

#include "statmap.h"
#include "mapped_read.h"
#include "find_candidate_mappings.h"
#include "index_genome.h"
#include "quality.h"
#include "iterative_mapping.h"
#include "candidate_mapping.h"
#include "fragment_length.h"
#include "sam.h"

/* fwd declaration */
struct snp_db_t;
#include "snp.h"

/* Set the deafults for these two global variables */
int num_threads = -1;
int min_num_hq_bps = -1;

void usage()
{
    fprintf(stderr, "Usage: ./statmap -g genome.fa -p men_match_penalty -m max_penalty_spread \n");
    fprintf(stderr, "                 ( -r input.fastq | [ -1 input.pair1.fastq & -2 input.pair2.fastq ] ) \n");
    fprintf(stderr, "      (optional)  [ -o output_directory  -a assay_type -n snp_cov -f fragment_lengths ] \n\n" );
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
        exit( -1 );
    }
}

input_file_type_t
guess_input_file_type( args_t* args )
{
    input_file_type_t input_file_type = UNKNOWN;

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
    for( i = 0; i < 1000; i++ )
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
guess_optimal_indexed_seq_len( args_t* args)
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
             "NOTICE      :  Setting Indexed Read Length to '%i' BPS\n", 
             seq_len );
    
    free_rawread( r );
    fclose( fp );
    
    return seq_len;
}

args_t
parse_arguments( int argc, char** argv )
{
    int error;

    /* 
     * initialize the structure that stores all of the 
     * configuration options.
     */
    args_t args;
    
    args.genome_fname = NULL;
    
    args.unpaired_reads_fnames = NULL;
    args.pair1_reads_fnames = NULL;
    args.pair2_reads_fnames = NULL;

    args.snpcov_fname = NULL;
    args.snpcov_fp = NULL;

    args.frag_len_fname = NULL;
    args.frag_len_fp = NULL;
    
    args.output_directory = NULL;
    
    args.sam_output_fname = NULL;
    
    args.log_fname = NULL;
    args.log_fp = NULL;
    
    args.min_match_penalty = 1;
    args.max_penalty_spread = -1;
    args.min_num_hq_bps = -1;

    args.num_starting_locations = 10;

    args.num_threads = -1;
    args.indexed_seq_len = -1;
    
    args.input_file_type = UNKNOWN;
    args.assay_type = UNKNOWN;

    char* assay_name = NULL;

    int c;
    while ((c = getopt(argc, argv, "9hg:n:r:1:2:c:o:p:m:s:f:l:a:w:t:q:")) != -1) {
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
            
        /* optional arguments ( that you probably want ) */
        case 'o': // output directory
            args.output_directory = optarg;
            break;

        case 'a': // the assay type
            assay_name = optarg;
            break;

        case 's': // snp input file
            args.snpcov_fname = optarg;
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
            
        /* optional arguments ( that you probably dont want ) */
        case 'i': // indexed sequence length
                  // defaults to the read length of the first read
            args.indexed_seq_len = atoi(optarg);
            break;
        case 'l': // log output 
            args.log_fname = optarg;
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
        fprintf(stderr, "FATAL       :  -g ( reference_genome ) is required\n");
        exit( -1 );
    }

    /* open the chromosome file */
    args.genome_fp = fopen( args.genome_fname, "r");
    if( args.genome_fp == NULL ) {
        fprintf(stderr, "FATAL       :  Unable to open '%s'\n", args.genome_fname);
        exit(-1);
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
        char buffer[500];
        sprintf( buffer, "cp %s %s/reads.unpaired", args.unpaired_reads_fnames, args.output_directory );
        fprintf(stderr, "NOTICE      :  Copying '%s' to the output directory\n",  args.unpaired_reads_fnames );
        error = system( buffer );
        if (WIFSIGNALED(error) &&
            (WTERMSIG(error) == SIGINT || WTERMSIG(error) == SIGQUIT))
        {
            fprintf(stderr, "FATAL     : Failed to call '%s'\n", buffer );
            perror( "System Call Failure");
            assert( false );
            exit( -1 );
        }
    } 
    /* If the reads are paired */
    else {
        /* first, copy the read file(s) into the output directory */
        char buffer[500];
        sprintf( buffer, "cp %s %s/reads.pair1", args.pair1_reads_fnames, args.output_directory );
        fprintf(stderr, "NOTICE      :  Copying '%s' to the output directory\n",  args.pair1_reads_fnames );
        error = system( buffer );
        if (WIFSIGNALED(error) &&
            (WTERMSIG(error) == SIGINT || WTERMSIG(error) == SIGQUIT))
        {
            fprintf(stderr, "FATAL     : Failed to call '%s'\n", buffer );
            perror( "System Call Failure");
            assert( false );
            exit( -1 );
        }

        fprintf(stderr, "NOTICE      :  Copying '%s' to the output directory\n",  args.pair2_reads_fnames );
        sprintf( buffer, "cp %s %s/reads.pair2", args.pair2_reads_fnames, args.output_directory );
        error = system( buffer );
        if (WIFSIGNALED(error) &&
            (WTERMSIG(error) == SIGINT || WTERMSIG(error) == SIGQUIT))
        {
            fprintf(stderr, "FATAL     : Failed to call '%s'\n", buffer );
            perror( "System Call Failure");
            assert( false );
            exit( -1 );
        }

    }

    /*
     * Try to determine the type of sequencing error. 
     */
    if( args.input_file_type == UNKNOWN )
    {
        args.input_file_type = guess_input_file_type( &args );
    }

    
    /*
     * If the sequence length is not set, then try and determine it automatically.
     */    
    if( args.indexed_seq_len == -1 )
    {
        args.indexed_seq_len = guess_optimal_indexed_seq_len( &args );
        
        if( args.indexed_seq_len <= 11 ) 
        {
            fprintf( stderr, "FATAL       :  Can not index sequences less than 12 basepairs long.\n" );
            exit( -1 );
        }
    }

    /* open the snp coverage file */
    if( args.snpcov_fname != NULL )
    {
        args.snpcov_fp = fopen( args.snpcov_fname, "r");
        if( NULL == args.snpcov_fp )
        {
            fprintf( stderr, "FATAL       :  Failed to open '%s'\n", args.snpcov_fname );
            exit( -1 );
        }
    }

    if( args.frag_len_fname != NULL )
    {
        args.frag_len_fp = fopen( args.frag_len_fname, "r" );
        if( NULL == args.frag_len_fp )
        {
            char buffer[500];
            sprintf( buffer, "FATAL       :  Failed to open '%s' ", args.frag_len_fname );
            perror( buffer );
            exit( -1 );
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

    if( SAVE_BOOTSTRAP_SAMPLES )
    {
        safe_mkdir( BOOTSTRAP_SAMPLES_PATH );
        safe_mkdir( BOOTSTRAP_SAMPLES_MAX_PATH );
        safe_mkdir( BOOTSTRAP_SAMPLES_MIN_PATH );
        safe_mkdir( BOOTSTRAP_SAMPLES_ALL_PATH );
    }

    /***** initialize the raw reads db */
    init_rawread_db( &(args.rdb) );

    /* If the reads are not paired */
    if( args.unpaired_reads_fnames != NULL )
    {
        add_single_end_reads_to_rawread_db(
            args.rdb, "reads.unpaired", FASTQ 
        );
    } 
    /* If the reads are paired */
    else {
        add_paired_end_reads_to_rawread_db(
            args.rdb, 
            "reads.pair1",
            "reads.pair2",
            FASTQ 
        );
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
        fprintf(stderr, "NOTICE      :  Number of threads is being set to %i \n", num_threads);
    } else {
        num_threads = args.num_threads;
    }


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

    /* initialize the log file */
    if( args.log_fname != NULL )
        args.log_fp = fopen( args.log_fname, "a");
    
    /* set the min num hq basepairs if it's unset */
    if( args.min_num_hq_bps == -1 )
    {
        args.min_num_hq_bps = MAX( 12, (8 + args.indexed_seq_len/2) );
    }
    min_num_hq_bps = args.min_num_hq_bps;

    /* make sure that we dont filter *every* read */
    if( min_num_hq_bps > args.indexed_seq_len )
    {
        fprintf( stderr, "WARNING     :  Required number of HQ bps set higher than the seq length. Setting it to seq length.\n" );
        min_num_hq_bps = args.indexed_seq_len;
    }

    /* 
     * Dont allow penalty spreads greater than the min match penalty - 
     * they are worthless and mess up my search heuristics 
     *
     */
    if( args.max_penalty_spread + args.min_match_penalty + 0.00001 > 0 )
        args.max_penalty_spread = -1;
    
    return args;
}


void
map_marginal( args_t* args, 
              struct genome_data* genome, 
              struct mapped_reads_db** mpd_rds_db )
{
    /* store clock times - useful for benchmarking */
    clock_t start, stop;
    
    /***** Index the genome */
    start = clock();
    
    /* initialize the genome */
    index_genome( genome, args->indexed_seq_len );
    
    stop = clock();
    fprintf(stderr, "PERFORMANCE :  Indexed Genome in %.2lf seconds\n", 
                    ((float)(stop-start))/CLOCKS_PER_SEC );
    
    fprintf(stderr, "PERFORMANCE :  Tree Size: %lu bytes\n", 
            (unsigned long) sizeof_tree( genome->index ) );
    
    /***** END Genome processing */
    
    
    /***** initialize the mappings dbs */
    
    candidate_mappings_db candidate_mappings;
    init_candidate_mappings_db( &candidate_mappings, "candidate_mappings" );
    
    init_mapped_reads_db( mpd_rds_db, MAPPED_READS_DB_FNAME );

    if( args->frag_len_fp != NULL )
        build_fl_dist_from_file( *mpd_rds_db, args->frag_len_fp );

    /***** END initialize the mappings dbs */
    
    start = clock();
    
    find_all_candidate_mappings( 
        genome,
        args->log_fp,
        args->rdb,
        &candidate_mappings,
        args->min_match_penalty,
        args->max_penalty_spread,
        args->indexed_seq_len );
    
    /* Free the genome index */
    /* we may need the memory later */
    fprintf(stderr, "NOTICE      :  Freeing index\n" );
    fprintf(stderr, "PERFORMANCE :  Tree Size: %lu bytes\n", 
            (unsigned long) sizeof_tree( genome->index ) );
    free_tree( genome->index );
    genome->index = NULL;
    
    /* combine and output all of the partial mappings - this includes
       joining paired end reads. */
    fprintf(stderr, "NOTICE      :  Joining Candidate Mappings\n" );
    start = clock();
    join_all_candidate_mappings( &candidate_mappings, *mpd_rds_db );
    stop = clock();
    fprintf(stderr, "PERFORMANCE :  Joined Candidate Mappings in %.2lf seconds\n", 
                    ((float)(stop-start))/CLOCKS_PER_SEC );
    /*  close candidate mappings db */
    close_candidate_mappings_db( &candidate_mappings );

    /* Write the mapped reads to file */
    fprintf(stderr, "NOTICE      :  Writing mapped reads to wiggle file.\n" );
    FILE* fwd_wig_fp = fopen( "marginal_mappings_fwd.wig", "w+" );
    FILE* bkwd_wig_fp = fopen( "marginal_mappings_bkwd.wig", "w+" );
    write_marginal_mapped_reads_to_stranded_wiggles( 
        *mpd_rds_db, genome, fwd_wig_fp, bkwd_wig_fp );
    fclose( fwd_wig_fp );
    fclose( bkwd_wig_fp );

    /* write the mapped reads to SAM */
    start = clock();
    fprintf(stderr, "NOTICE      :  Writing non mapping reads to FASTQ files.\n" );
    write_nonmapping_reads_to_fastq( args->rdb, *mpd_rds_db );
    stop = clock();
    fprintf(stderr, "PERFORMANCE :  Wrote non-mapping reads to FASTQ in %.2lf sec\n", 
                    ((float)(stop-start))/CLOCKS_PER_SEC );

    /* write the mapped reads to SAM */
    start = clock();
    fprintf(stderr, "NOTICE      :  Writing mapped reads to SAM file.\n" );
    FILE* sam_ofp = fopen( "mapped_reads.sam", "w+" );
    write_mapped_reads_to_sam( 
        args->rdb, *mpd_rds_db, genome, false, false, sam_ofp );
    fclose( sam_ofp );    
    stop = clock();
    fprintf(stderr, "PERFORMANCE :  Wrote mapped reads to sam in %.2lf seconds\n", 
                    ((float)(stop-start))/CLOCKS_PER_SEC );
    
    return;
}

void
build_fl_dist( args_t* args, struct mapped_reads_db* mpd_rds_db )
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
iterative_mapping( args_t* args, 
                   struct genome_data* genome, 
                   struct mapped_reads_db* mpd_rds_db )
{   
    if( NULL != args->unpaired_reads_fnames 
        && args->assay_type == CHIP_SEQ )
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
    generic_update_mapping( args->rdb, 
                            mpd_rds_db,
                            genome, 
                            args->assay_type,
                            args->num_starting_locations, 
                            MAX_PRB_CHANGE_FOR_CONVERGENCE );
        
    munmap_mapped_reads_db( mpd_rds_db );
    
    return;
}

int 
main( int argc, char** argv )
{       
    /* parse and sanity check arguments */
    args_t args = parse_arguments( argc, argv );
    
    /* Load the genome */
    struct genome_data* genome;
    init_genome( &genome );
    add_chrs_from_fasta_file( genome, args.genome_fp );

    /* parse the snps */
    if( args.snpcov_fp != NULL )
        build_snp_db_from_snpcov_file( args.snpcov_fp, genome );

    struct mapped_reads_db* mpd_rds_db;
           
    /* Map the marginal reads and output them into a read density wiggle */
    map_marginal( &args, genome, &mpd_rds_db );

    build_fl_dist( &args, mpd_rds_db );

    if( args.assay_type != UNKNOWN )
        iterative_mapping( &args, genome, mpd_rds_db );

    /* If appropriate, print out the snp db */
    if( args.snpcov_fp != NULL )
    {
        fprintf(stderr, "NOTICE      :  Updating SNP count estiamtes.\n" );
        update_snp_estimates_from_candidate_mappings( mpd_rds_db, genome );
        char* snp_fname = "updated_snp_cnts.snp";
        FILE* snp_fp = fopen( snp_fname, "w" );
        write_snps_to_file( snp_fp, genome );
        fclose( snp_fp );
    }
    
    goto cleanup;
    
cleanup:
    /* Free the genome and indexes */
    free_genome( genome );
    
    /* Close the log file */
    if( args.log_fp != NULL )
        fclose(args.log_fp);
    
    /* close the snp coverage file */
    if( args.snpcov_fp != NULL ) {
        fclose( args.snpcov_fp );
    }

    return 0;
}
