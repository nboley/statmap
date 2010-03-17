/* Copyright (c) 2009-2010, Nathan Boley */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <getopt.h>
#include <string.h>
#include <pthread.h>

#include "statmap.h"
#include "find_candidate_mappings.h"
#include "index_genome.h"
#include "db_interface.h"
#include "quality.h"
#include "iterative_mapping.h"

/* fwd declaration */
struct snp_db_t;
#include "snp.h"

/* Store parsed options */
typedef struct {
    char* genome_fname;

    char* unpaired_reads_fnames;
    char* pair1_reads_fnames;
    char* pair2_reads_fnames;

    char* snpcov_fname;
    FILE* snpcov_fp;

    char* wig_fname;
    FILE* wig_fp;
    
    char* candidate_mappings_prefix;

    char* sam_output_fname;

    char* log_fname;
    FILE* log_fp;

    float min_match_penalty;
    float max_penalty_spread;

    int indexed_seq_len;
    
    input_file_type_t input_file_type;
    enum assay_type_t assay_type; 
} args_t;

void usage() 
{
    fprintf(stderr, "Usage: ./statmap -g genome.fa -p men_match_penalty -m max_penalty_spread \n");
    fprintf(stderr, "                 ( -r input.fastq | [ -1 input.pair1.fastq & -2 input.pair2.fastq ] ) \n");
    fprintf(stderr, "      (optional) [ -o output.sam  -s indexed_seq_length -l logfile ] \n\n" );
}

/*
 * Try and determine the file type. 
 *
 * In particular, we try and determine what the sequence mutation
 * string types are.
 *
 * The method is to scan the first 10000 reads and record the 
 * max and min untranslated scores. 
 *
 *
 *
 */

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
        rawread* r;
        populate_read_from_fastq_file( fp, &r );
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
        if( max_qual == 104 ) {
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
            else if ( min_qual < 69 ) {
                fprintf(stderr, "WARNING     :  input file format ambiguous.\n");
                input_file_type = SOLEXA_LOG_ODDS_FQ;
            } else {
                fprintf( stderr, "ERROR       :  Could not automatically determine " );
                fprintf( stderr, "input format. ( max=%i, min=%i  )\n", 
                         max_qual, min_qual );
                fprintf( stderr, "ERROR       :  SETTING TO %i. Results could be incorrect.\n", ILLUMINA_v13_FQ );
                input_file_type = ILLUMINA_v13_FQ;
                // exit( -1 );
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

/*
 * Guess the optimal indexed sequence length.
 *
 * To do this, we open up the first read file and then scan for a read. The 
 * seq length of the first read is what we set the index length to.
 *
 */
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
    rawread* r;
    populate_read_from_fastq_file( fp, &r );
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

    args.wig_fname = NULL;
    args.wig_fname = NULL;
    
    args.candidate_mappings_prefix = NULL;
    args.sam_output_fname = NULL;

    args.log_fname = NULL;
    args.log_fp = NULL;
    
    args.min_match_penalty = -1;
    args.max_penalty_spread = -1;
    
    args.indexed_seq_len = -1;
    
    args.input_file_type = UNKNOWN;
    args.assay_type = UNKNOWN;

    char* assay_name = NULL;

    int c;
    while ((c = getopt(argc, argv, "9hg:n:r:1:2:c:o:p:m:s:l:a:w:")) != -1) {
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
        case 'o': // output file name 
            args.sam_output_fname = optarg;
            break;
        case 'a': // the assay type
            assay_name = optarg;
            break;
        case 'n': // snp input file
            args.snpcov_fname = optarg;
            break;
        case 'w': // wig ouput file
            args.wig_fname = optarg;
            break;

        /* optional arguments ( that you probably dont want ) */
        case 's': // indexed sequence length
                  // defaults to the read length of the first read
            args.indexed_seq_len = atoi(optarg);
            break;
        case 'c': // directory to store candidate mappings in 
                  // during the mapping process, statmap stores
                  // all of the partial mappings in a list of 
                  // files - defaults to a random directory
            args.candidate_mappings_prefix = optarg;
            break;
        case 'l': // log output 
            args.log_fname = optarg;
            break;

        /* utility options */
        case 'h':
            usage();
            exit(-1);            
        case '?':
            fprintf(stderr, "ERROR       :  Unrecognized Argument: '%c' \n", (char) optopt);
            usage();
            exit(-1);
        default:
            usage();
            exit( -1 );
        }
    }

    /* Ensure that the required arguments were present */
    if( args.genome_fname == NULL )
    {
        usage();
        fprintf(stderr, "ERROR       :  -g ( reference_genome ) is required\n");
        exit( -1 );
    }

    if( args.min_match_penalty == -1 )
    {
        usage();
        fprintf(stderr, "ERROR       :  -p ( min match penalty ) is required\n");
        exit( -1 );
    }

    if( args.max_penalty_spread == -1 )
    {
        usage();
        fprintf(stderr, "ERROR       :  -m ( max penalty spread ) is required\n");
        exit( -1 );
    }

    if( args.unpaired_reads_fnames == NULL
        && ( args.pair2_reads_fnames == NULL
             || args.pair1_reads_fnames == NULL ) 
        )
    {
        usage();
        fprintf(stderr, "ERROR       :  -r or ( -1 and -2 ) is required\n" );
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
             "ERROR       :  if the -r is set, neither -1 nor -2 should be set\n"
            );
            exit(-1);
        }
    }

    if( args.pair1_reads_fnames != NULL 
        && args.pair2_reads_fnames == NULL )
    {
        usage();
        fprintf(stderr, 
                "ERROR       :  -1 option is set but -2 is not\n" );
        exit(-1);            
    }

    if( args.pair2_reads_fnames != NULL 
        && args.pair1_reads_fnames == NULL )
    {
        usage();
        fprintf(stderr, 
                "ERROR       :  -2 option is set but -1 is not\n" );
        exit(-1);            
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

    /* open the snp coverage file */
    if( args.snpcov_fname != NULL )
    {
        args.snpcov_fp = fopen( args.snpcov_fname, "r");
        if( NULL == args.snpcov_fp )
        {
            fprintf( stderr, "FATAL       :  Failed to open '%s'\n", args.snpcov_fname );
            exit( 1 );
        }
    }

    /* open the wig file */
    if( args.wig_fname != NULL )
    {
        args.wig_fp = fopen( args.wig_fname, "w");
        if( NULL == args.wig_fp )
        {
            fprintf( stderr, "FATAL       :  Failed to open '%s' for writing\n", args.wig_fname );
            exit( 1 );
        }
    }
    
    /*
     * If the sequence length is not set, then try and determine it automatically.
     */
    
    if( args.indexed_seq_len == -1 )
    {
        args.indexed_seq_len = guess_optimal_indexed_seq_len( &args );
    }

    /*
     * Try to determine the type of sequencing error. 
     */
    if( args.input_file_type == UNKNOWN )
    {
        args.input_file_type = guess_input_file_type( &args );
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
join_all_candidate_mappings( candidate_mappings_db* cand_mappings_db,
                             mapped_reads_db* mpd_rds_db )
{
    int error;

    unsigned int read_id = 0;
    
    clock_t start, stop;
    start = clock();
    
    /* Join all candidate mappings */

    /* get the cursor to iterate through the candidate mappings */    
    candidate_mappings_db_cursor* candidate_mappings_cursor;
    open_candidate_mappings_cursor(
        cand_mappings_db, &candidate_mappings_cursor );
    
    /* BUG - what is this here for ? */
    char curr_key[ MAX_KEY_SIZE + 1];
    
    candidate_mappings* mappings;
    mapped_read* mpd_rd;

    /* Get the first read */
    error = get_next_candidate_mapping_from_cursor( 
        candidate_mappings_cursor, 
        &mappings,
        curr_key 
    );
    
    while( CURSOR_EMPTY != error ) 
    {
        build_mapped_read_from_candidate_mappings( 
            mappings, &mpd_rd, read_id );
        
        add_read_to_mapped_reads_db( mpd_rds_db, mpd_rd );

        free_mapped_read( mpd_rd );

        free_candidate_mappings( mappings );

        /* Get the reads */
        error = get_next_candidate_mapping_from_cursor( 
            candidate_mappings_cursor, 
            &mappings,
            curr_key 
        );

        read_id += 1;
    }
    
    goto cleanup;
    
cleanup:
    /* close the cursor */
    close_candidate_mappings_cursor( candidate_mappings_cursor );
    
    stop = clock();
    fprintf(stderr, "PERFORMANCE :  Joined Candidate Mappings in %.2lf seconds\n", 
                    ((float)(stop-start))/CLOCKS_PER_SEC );
    
    return;
}


void
write_mapped_reads_to_sam( rawread_db_t* rdb,
                                  mapped_reads_db* mappings_db,
                                  genome_data* genome,
                                  FILE* sam_ofp )
{
    int error;
    
    clock_t start, stop;
    start = clock();
    
    /* Join all candidate mappings */
    /* get the cursor to iterate through the reads */
    rewind_rawread_db( rdb );
    rewind_mapped_reads_db( mappings_db );
    
    rawread *rd1, *rd2;
    mapped_read* mapped_rd;

    error = get_next_read_from_mapped_reads_db( 
        mappings_db, 
        &mapped_rd
    );
    
    get_next_mappable_read_from_rawread_db( 
        rdb, &rd1, &rd2 );
    
    while( !mapped_reads_db_is_empty( mappings_db ) ) 
    {     
        /* 
         * If we couldnt map it anywhere,
         * print out the read 
         */
        if( mapped_rd->num_mappings == 0 )
        {
            /* If this is a single end read */
            if( rd2 == NULL )
            {
                fprintf_rawread_to_fastq( 
                    rdb->non_mapping_single_end_reads, rd1 );
            } else {
                fprintf_rawread_to_fastq( 
                    rdb->non_mapping_paired_end_1_reads, rd1 );

                fprintf_rawread_to_fastq( 
                    rdb->non_mapping_paired_end_2_reads, rd2 );
            }
        /* otherwise, print it out to the sam file */
        } else {
            fprintf_mapped_read_to_sam( 
                sam_ofp, mapped_rd, genome, rd1, rd2 );
        }

        free_mapped_read( mapped_rd );
        
        /* Free the raw reads */
        free_rawread( rd1 );
        if( rd2 != NULL )
            free_rawread( rd2 );

        error = get_next_read_from_mapped_reads_db( 
            mappings_db, 
            &mapped_rd
        );

        get_next_mappable_read_from_rawread_db( 
            rdb, &rd1, &rd2 );
        
    }

    goto cleanup;

cleanup:
    free_mapped_read( mapped_rd );
        
    stop = clock();
    fprintf(stderr, "PERFORMANCE :  Wrote mapped reads to sam in %.2lf seconds\n", 
                    ((float)(stop-start))/CLOCKS_PER_SEC );
    
    return;
}

void
map_marginal( args_t* args, genome_data* genome )
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
    
    candidate_mappings_db mappings_db;
    init_candidate_mappings_db( &mappings_db, 
                                args->candidate_mappings_prefix );
    
    mapped_reads_db* mpd_rds_db;
    init_mapped_reads_db( &mpd_rds_db, "test.mapped_reads_db" );

    /***** END initialize the mappings dbs */
    
    /***** initialize the raw reads db */

    rawread_db_t* raw_rdb;
    init_rawread_db( &raw_rdb );

    /* If the reads are not paired */
    if( args->unpaired_reads_fnames != NULL )
    {
        add_single_end_reads_to_rawread_db(
            raw_rdb, args->unpaired_reads_fnames, FASTQ 
        );
    } 
    /* If the reads are paired */
    else {
        add_paired_end_reads_to_rawread_db(
            raw_rdb, 
            args->pair1_reads_fnames, 
            args->pair2_reads_fnames, 
            FASTQ 
        );
    }
    /***** END initialize the read db */

    start = clock();
    
    find_all_candidate_mappings( genome,
                                 args->log_fp,
                                 raw_rdb,
                                 &mappings_db,
                                 args->min_match_penalty,
                                 args->max_penalty_spread,
                                 args->indexed_seq_len );
    
    /* Determine the output stream */
    FILE* sam_ofp = stdout;
    if( args->sam_output_fname != NULL )
        sam_ofp = fopen( args->sam_output_fname, "w" );        
    
    /* combine and output all of the partial mappings - this includes
       joining paired end reads. */
    join_all_candidate_mappings( &mappings_db, mpd_rds_db );
    
    /* Iterative mapping */
    /* mmap and index the necessary data */
    mmap_mapped_reads_db( mpd_rds_db );
    index_mapped_reads_db( mpd_rds_db );
    update_mapping( mpd_rds_db, genome, 50, args->assay_type  );
        
    if( args->wig_fp != NULL )
        write_mapped_reads_to_wiggle( mpd_rds_db, genome, args->wig_fp );

    munmap_mapped_reads_db( mpd_rds_db );
    
    write_mapped_reads_to_sam( 
        raw_rdb, mpd_rds_db, genome, sam_ofp );

    goto cleanup;

cleanup:
    /* Free the genome index */
    free_tree( genome->index );
    genome->index = NULL;
        
    /* close the raw mappings db */
    close_rawread_db( raw_rdb );
    
    /*  close candidate mappings db */
    /*  This removes the temporary files as well */
    close_candidate_mappings_db( &mappings_db );

    /* Close the packed mapped reads db */
    close_mapped_reads_db( mpd_rds_db );

    /* if we opened an output file, close it */
    if( args->sam_output_fname != NULL )
        fclose( sam_ofp );    
    
    return;
}

int 
main( int argc, char** argv )
{          
    /* parse and sanity check arguments */
    args_t args = parse_arguments( argc, argv );
    
    /* Load the genome */
    genome_data* genome;
    init_genome( &genome );
    add_chrs_from_fasta_file( genome, args.genome_fname );

    /* parse the snps */
    if( args.snpcov_fp != NULL )
        build_snp_db_from_snpcov_file( args.snpcov_fp, genome );
           
    /* Map the marginal reads and output them into a sam */
    map_marginal( &args, genome );
    
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

    /* close the wig file */
    if( args.wig_fp != NULL ) {
        fclose( args.wig_fp );
    }

    return 0;
}
