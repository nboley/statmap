/* Copyright (c) 2009-2012, Nathan Boley */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <assert.h>

#include <time.h>
#include <sys/time.h> /* gettimeofday() */

/* to find out the  number of available processors */
#include <sys/sysinfo.h>
/* mkdir */
#include <sys/stat.h>
/* chdir */
#include <sys/unistd.h>
/* check for system signals */
#include <signal.h>
#include <sys/wait.h>

#include <argp.h>

#include "config_parsing.h"
#include "log.h"
#include "find_candidate_mappings.h"
#include "rawread.h"
#include "quality.h"
#include "statmap.h"
#include "error_correction.h"
#include "util.h"

/**** Handle different sources of reads ****/

void
set_global_quality_parameters_from_input_file_type(
        enum input_file_type_t input_file_type
    )
{
    /* print out the determined format */
    statmap_log( LOG_NOTICE, "Setting input file format to %i",
            input_file_type );

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
#if 0
    default:
        statmap_log( LOG_FATAL, "Unrecognized file format type %i",
                input_file_type );
#endif
    }

    statmap_log( LOG_NOTICE, "Set QUAL_SHIFT to '%i' and ARE_LOG_ODDS to %i",
            QUAL_SHIFT, ARE_LOG_ODDS );
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
        populate_rawread_from_fastq_file( fp, &r, NORMAL );
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
    fclose( fp );

    statmap_log( LOG_NOTICE, "Calculated max and min quality scores as '%i' and '%i'",  max_qual, min_qual  );

    /* if the max quality score is 73, ( 'I' ), then
       this is almost certainly a sanger format */
    if( max_qual == 73 ) {
        input_file_type = SANGER_FQ ;
    } else if( max_qual == 74 ) {
        /* if the max quality score is 74, ('J'), then
         * this is almost certainly Illumina 1.8+ */
        input_file_type = ILLUMINA_v18_FQ;
    } else if( max_qual > 74 && max_qual <= 90 ) {
        input_file_type = MARKS_SOLEXA;
    } else {
        // If the max quality is 104 it's a bit tougher.
        // because it can be either the new or old illumina format
	if( max_qual < 104 )
	{
	    /* If the maximum quality is less than 104 but greater than 73, it's probably the
	       plus 64 version, but we can't be sure. However, assume that it is and print a 
	       warning. 
	    */
	    statmap_log( LOG_WARNING, "Maximum input score wasn't achieved. ( %i )", max_qual);
	    statmap_log( LOG_WARNING, "( That means the highest quality basepair was above 73 but less than 104. Probably, there are just no very HQ basepairs, but make sure that the predicted error format is correct. )" );
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
	    statmap_log( LOG_WARNING, "input file format ambiguous. ( max=%i, min=%i )",  max_qual, min_qual );
	    input_file_type = ILLUMINA_v13_FQ;
	} else if ( min_qual < 70 ) {
	    statmap_log( LOG_WARNING, "input file format ambiguous. ( max=%i, min=%i )",  max_qual, min_qual );
	    input_file_type = ILLUMINA_v15_FQ;
	} else if ( min_qual < 90 ) {
	    statmap_log( LOG_WARNING, "input file format ambiguous. ( max=%i, min=%i )",  max_qual, min_qual );
	    input_file_type = MARKS_SOLEXA;
	} else if ( min_qual == 104 ) {
	    statmap_log( LOG_WARNING, "input file format ambiguous. ( max=%i, min=%i )",  max_qual, min_qual );
	    input_file_type = TEST_SUITE_FORMAT;
	} else {
	    statmap_log( LOG_FATAL, "Could not automatically determine input format. ( max=%i, min=%i  )",  max_qual, min_qual  );
	}
    }

    set_global_quality_parameters_from_input_file_type( input_file_type );
    return input_file_type;
}

int
set_input_file_type( char* arg, struct args_t* args )
{
    int flag = atoi( arg );

    switch( flag )
    {
        case 1:
            args->input_file_type = SANGER_FQ;
            break;
        case 2:
            args->input_file_type = ILLUMINA_v13_FQ;
            break;
        case 3:
            args->input_file_type = ILLUMINA_v15_FQ;
            break;
        case 4:
            args->input_file_type = SOLEXA_v14_FQ;
            break;
        case 5:
            args->input_file_type = SOLEXA_LOG_ODDS_FQ;
            break;
        case 6:
            args->input_file_type = TEST_SUITE_FORMAT;
            break;
        case 7:
            args->input_file_type = MARKS_SOLEXA;
            break;
        case 8:
            args->input_file_type = ILLUMINA_v18_FQ;
            break;
        default:
            /* the flag type was unrecognized; signal failure with -1 */
            return -1;
    }

    /* if the flag was recognized, set the appropriate global variables */
    set_global_quality_parameters_from_input_file_type( args->input_file_type );

    /* return the value of the successfully passed input file type flag */
    return flag;
}

/******** Parse command line arguments ********/

const char* argp_program_version = 
"Statmap 0.3.3";
const char* argp_program_bug_address =
"npboley@gmail.com";

/*
   OPTIONS.  Field 1 in ARGP.
   Order of fields: {NAME, KEY, ARG, FLAGS, DOC, GROUP}.

   Note - if entry has GROUP = 0, inherits previous entry's GROUP
   If GROUP is not specified, it defaults to 0.

   Can add options purely for documentation purposes if all options except DOC
   are 0. These entries auto-increment the GROUP.
*/
static struct argp_option options[] = 
{
    /* genome file -- required */
    {0, 0, 0, 0, "Genome:", 0},
    {"genome", 'g', "GENOME", 0,
     "Path to Statmap format genome file", 0},

    /* input reads -- required */
    {0, 0, 0, 0, "Single End Reads:", 0},
    {"unpaired", 'r', "READS", 0, 
     "FASTQ input file for single end reads", 0},
    {"unpaired-negative-control", 'c', "NCONTROL", 0,
     "Specifies negative control data (ChIP-Seq only, optional)", 0},

    {0, 0, 0, 0, "Paired End Reads:", 0},
    {"paired-1", '1', "READS1", 0,
     "FASTQ input for for the first set of read pairs", 0},
    {"paired-2", '2', "READS2", 0,
     "FASTQ input for for the second set of read pairs", 0},
    {"negative-control-1", '3', "NCONTROL1", 0,
     "Specifies pair 1 negative control data (ChIP-Seq only, optional)", 0},
    {"negative-control-2", '4', "NCONTROL2", 0,
     "Specifies pair 1 negative control data (ChIP-Seq only, optional)", 0},

    /* optional arguments */
    {0, 0, 0, 0, "Optional Arguments:", 0},
    {"mapping-metaparameter", 'p', "PARAM", 0,
     "Fraction of reads to try to map (estimated error model) or number of mismatches allowed as fraction of read length (mismatch error model)", 0},
    {"max-penalty-spread", 'm', "SPREAD", 0,
     "Upper bound of the difference of probabilities between the top matching sequence and the sequence of interest (in log10 probaiblity space)", 0},
    {"assay", 'a', "ASSAY", 0,
     "Optional: type of underlying assay. Valid options are 'i' for ChIP-Seq or 'a' for CAGE. 'r' for stranded RNA-Seq is not implemented yet", 0},
    {"num-samples", 'n', "SAMPLES", 0,
     "In iterative mapping, number of samples to take from the mapping posterior", 0},
    {"threads", 't', "THREADS", 0,
     "Number of threads to use. Defaults to all available, but no more than 8.", 0},
    {"frag-len-dist", 'f', "DIST", 0,
     "Fragment length distribution file", 0},
    {"output-dir", 'o', "DIR", 0,
     "Directory to write Statmap output to", 0 },
    {"search-type", 's', "TYPE", 0,
     "Error model type: 'm' for mismatches, 'e' to estimate the error model ( see README for more details )", 0 },
    {"is-full-fragment", 'F', NULL, OPTION_ARG_OPTIONAL,
     "(not implemented yet)", 0 },
    {"input-file-type", 'i', "TYPE", 0,
     "Type of sequencing platform used to generate the input reads. Valid "
     "options are 1: SANGER_FQ, 2: ILLUMINA_v13_FQ, 3: ILLUMINA_v15_FQ, 4: "
     "SOLEXA_v14_FQ, 5: SOLEXA_LOG_ODDS_FQ, 6: TEST_SUITE_FORMAT, 7: "
     "MARKS_SOLEXA, 8: ILLUMINA_v18_FQ", 0 },
    {"paired-end-reads-map-to-opposite-strands", 'P', NULL, 0,
     "(not implemented yet)", 0 },
    {"soft-clip-length", 'S', "LEN", 0,
     "Number of basepair's to soft clip from the start of each read (default 0)", 0 },
    {"verbose", 'V', NULL, OPTION_ARG_OPTIONAL,
     "Set the logging level to NOTICE", 0 },
    {"debug-verbose", 'D', NULL, OPTION_ARG_OPTIONAL,
     "Set the logging level to DEBUG", 0 },

    /* must end with an entry containing all zeros */
    {0,0,0,0,0,0}
};

/*
   PARSER. Field 2 in ARGP.
   Order of parameters: KEY, ARG, STATE.
*/
static error_t
parse_opt( int key, char *arg, struct argp_state *state )
{
    struct args_t* args = state->input;

    switch(key)
    {
        /* Required arguments - genome and reads */
        case 'g':
            args->genome_fname = arg;
            break;
        case 'r':
            args->unpaired_reads_fnames = arg;
            break;
        case '1':
            args->pair1_reads_fnames = arg;
            break;
        case '2':
            args->pair2_reads_fnames = arg;
            break;
        case 'c':
            args->unpaired_NC_reads_fnames = arg;
            break;
        case '3':
            args->pair1_NC_reads_fnames = arg;
            break;
        case '4':
            args->pair2_NC_reads_fnames = arg;
            break;

            /* Optional arguments */
        case 'p':
            args->mapping_metaparameter = atof(arg);
            break;
        case 'o':
            args->output_directory = arg;
            break;
        case 'a':
            /* set args->assay_type depending on argument */
            switch( arg[0] )
            {
                case 'a':
                    args->assay_type = CAGE;
                    break;
                case 'i':
                    args->assay_type = CHIP_SEQ;
                    break;
                case 'r':
                    args->assay_type = RNA_SEQ;
                    break;
                default:
                    argp_failure( state, 1, 0, 
                            "FATAL       :  Unrecognized assay type: '%s'",
                            arg );
            }

            /* Set the global assay type variable */
            _assay_type = args->assay_type;

            break;
        case 'f':
            args->frag_len_fname = arg;
            break;
        case 'n':
            args->num_starting_locations = atoi(arg);
            break;
        case 't':
            args->num_threads = atoi(arg);
            break;
        case 's':
            switch( arg[0] )
            {
                case 'e':
                    args->error_model_type = ESTIMATED;
                    break;
                case 'm':
                    args->error_model_type = MISMATCH;
                    break;
                default:
                    argp_failure( state, 1, 0,
                            "FATAL       :  Invalid search type '%s'",
                            arg );
            }
            break;
        case 'F':
            argp_failure( state, 1, 0,
                    "FATAL       :  -F ( --is-full-fragment ) is not implemented yet." );
            break;
        case 'i':
        {
            // Note - cannot declare variables inside of a case without {}
            int rv = set_input_file_type( arg, args );
            if( rv < 0 )
            {
                argp_failure( state, 1, 0,
                        "FATAL        :  Unrecognized input file type '%s'",
                        arg );
            }
            break;
        }
        case 'P':
            argp_failure( state, 1, 0,
                    "FATAL       :  -P ( --paired-end-reads-map-to-opposite-strands ) is not implemented yet." );
            break;
        case 'S':
            args->softclip_len = atoi( arg );
            break;
        case 'V':
            args->nontrivial_log_level = LOG_INFO;
            break;
        case 'D':
            args->nontrivial_log_level = LOG_DEBUG;
            break;

            /* utility options */
            /* --help and --version are automatically provied by argp */
        default:
            return ARGP_ERR_UNKNOWN;
    }

    return 0;
}

/*
   ARGS_DOC. Field 3 in ARGP.
   A description of the non-option command-line arguments
*/
static char args_doc[] = "";

/*
   DOC.  Field 4 in ARGP.
   Program documentation.
*/
static char doc[] =
"Statmap is a fast alignment tool for mapping short reads to a reference.";

/* the ARGP structure */
static struct argp argp = {options, parse_opt, args_doc, doc, 0, 0, 0};

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

    args.mapping_metaparameter = -1;

    args.num_starting_locations = -1;

    args.num_threads = -1;

    args.error_model_type = UNKNOWN;
    args.input_file_type = UNKNOWN;
    args.assay_type = UNKNOWN;

    args.max_reference_insert_len = -1;
    args.softclip_len = 0;

    args.nontrivial_log_level = LOG_NOTICE;

    /* parse arguments with argp */
    argp_parse( &argp, argc, argv, 0, 0, &args );

    /* check config, open fps, etc. */

    /********* CHECK REQUIRED ARGUMENTS *************************************/

    if( args.genome_fname == NULL ) {
        statmap_log( LOG_FATAL, "-g (binary genome) is required" );
        exit(-1);
    }

    if( args.unpaired_reads_fnames == NULL
            && ( args.pair1_reads_fnames == NULL ||
                args.pair2_reads_fnames == NULL ))
    {
        statmap_log( LOG_FATAL, "-r or (-1 and -2) is required" );
        exit(-1);
    }

    /* perform sanity checks on the read input arguments */
    if( args.unpaired_reads_fnames )
    {
        if( args.pair1_reads_fnames != NULL || args.pair2_reads_fnames != NULL )
        {
            statmap_log( LOG_FATAL, "if -r is set, neither -1 nor -2 should be set" );
        }
    }

    if( args.pair1_reads_fnames != NULL &&
            args.pair2_reads_fnames == NULL )
    {
        statmap_log( LOG_FATAL, "-1 is set but -2 is not" );
    }

    if( args.pair2_reads_fnames != NULL &&
            args.pair1_reads_fnames == NULL )
    {
        statmap_log( LOG_FATAL, "-2 is set but -1 is not" );
    }

    /* If no output directory specified, make it */

    /* Build the default output directory filename from the current time */
    if( args.output_directory == NULL )
    {
        time_t rawtime;
        time( &rawtime );
        struct tm* timeinfo;
        timeinfo = localtime( &rawtime );

        char buffer[200];
        strftime( buffer, 200, "statmap_output_%Y_%m_%d_%H_%M_%S", timeinfo );
        args.output_directory = calloc( strlen(buffer)+1, sizeof(char) );
        strncpy( args.output_directory, buffer, strlen(buffer) );
    }
    /* copy the output directory argument into a new buffer so output_directory
     * can be freed later */
    else {
        char* saved_ptr = args.output_directory; // this is a char* into argv
        args.output_directory = calloc( strlen(args.output_directory)+1,
                sizeof(char) );
        strncpy( args.output_directory, saved_ptr, strlen(saved_ptr) );
    }

    error = mkdir( args.output_directory, S_IRWXU | S_IRWXG | S_IRWXO );
    if( -1 == error )
    {
        statmap_log( LOG_FATAL, "Cannot make output directory" );
        exit( -1 );
    }

    /* Check path of genome file and expand to absolute paths */
    char* genome_fname = realpath( args.genome_fname, NULL );
    if( NULL == genome_fname )
    {
        statmap_log( LOG_FATAL, "Could not determine absolute path of genome file" );
        exit( -1 );
    }
    args.genome_fname = genome_fname;

    /* Set args.genome_index_fname based on args.genome_fname */
    args.genome_index_fname = calloc( PATH_MAX - 6, sizeof(char) );
    sprintf( args.genome_index_fname, "%s.index", args.genome_fname );

    /* Set default for search type (necessary so we know what defaults to use
     * for mapping parameters) */
    if( args.error_model_type == UNKNOWN )
    {
        args.error_model_type = ESTIMATED;
    }

    /* Set defaults for numeric parameters */

    if( args.mapping_metaparameter == -1 )
    {
        if( args.error_model_type == ESTIMATED )
        {
            args.mapping_metaparameter = DEFAULT_ESTIMATED_ERROR_METAPARAMETER;
            statmap_log( LOG_NOTICE, "Setting the mapping metaparameter (-p) (for the estimated error model) to %.3f",  DEFAULT_ESTIMATED_ERROR_METAPARAMETER  );
        }
        else if( args.error_model_type == MISMATCH )
        {
            args.mapping_metaparameter = DEFAULT_MISMATCH_METAPARAMTER;
            statmap_log( LOG_NOTICE, "Setting the mapping metaparameter (-p) (for the mismatch error model) to %.3f",  DEFAULT_MISMATCH_METAPARAMTER  );
        }
    }

    /* Make sure the metaparameter value is valid. We expect a fraction in the
     * range [0,1] */
    if( args.mapping_metaparameter< 0 || args.mapping_metaparameter > 1 )
    {
        statmap_log( LOG_FATAL, "The mapping metaparameter (-p) must be in the range [0,1] got %.3f",  args.mapping_metaparameter  );
        assert( false );
        exit(1);
    }

    /***** Copy the reads file into the output directory ****/

    /* If the reads are not paired */
    if( args.unpaired_reads_fnames != NULL )
    {
        /* first, copy the read file(s) into the output directory */
        safe_copy_into_output_directory( 
                args.unpaired_reads_fnames, args.output_directory,
                "reads.unpaired" );

        /* copy the unpaired NC reads */
        if( args.unpaired_NC_reads_fnames != NULL )
        {
            safe_copy_into_output_directory( 
                    args.unpaired_NC_reads_fnames, args.output_directory,
                    "reads.NC.unpaired" );
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
                    args.pair1_NC_reads_fnames, args.output_directory,
                    "reads.NC.pair1" );

            if( args.pair2_NC_reads_fnames == NULL )
            {
                statmap_log( LOG_FATAL, "The NC reads must be paired for a paired experiment." );
            }

            safe_copy_into_output_directory( 
                    args.pair2_NC_reads_fnames, args.output_directory,
                    "reads.NC.pair2" );
        }
    }

    /*
       Try to determine the type of sequencing error. 
       */
    if( args.input_file_type == UNKNOWN )
    {
        args.input_file_type = guess_input_file_type( &args );
    }

    /*
       Copy fragment length distribution
       If this a  ChIP-Seq assay, error if no fl dist provided 
    */
    if( args.frag_len_fname != NULL )
    {
        safe_copy_into_output_directory( 
                args.frag_len_fname, args.output_directory, "estimated_fl_dist.txt" );

        args.frag_len_fp = fopen( args.frag_len_fname, "r" );
        if( NULL == args.frag_len_fp )
        {
            statmap_log( LOG_FATAL, "Failed to open '%s'\n",  args.frag_len_fname );
            exit( -1 );
        }    
    } else {
        if( args.assay_type == CHIP_SEQ )
        {
            statmap_log( LOG_FATAL, "Chip-Seq assay mapping requires a fragment lenght distribution (-f)" );
        }
    }

    /********* END CHECK REQUIRED ARGUMENTS ***********************************/

    /* change the working directory to the output directory */
    error = chdir( args.output_directory );
    if( -1 == error )
    {
        statmap_log( LOG_FATAL, "Cannot move into output directory" );
        exit( -1 );
    }

    /* Now that we've chdir'ed into the output directory , we can initialize
     * the real logfile. */
    init_logging( args.nontrivial_log_level );

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
    /***** END initialize the raw reads db */

    /* If we didnt set the number of threads, set it to the maximum available */
    if( args.num_threads == -1 )
    {
        /* try to get the number of available threads from the os */
        args.num_threads = get_nprocs();
        /* if we cant determine the number of threads, set it to 1 */
        if( args.num_threads <= 0 )
            args.num_threads = 1;
        
        /* never set the number of threads to more than 8, by default */
        if( args.num_threads > 8 )
            args.num_threads = 8;
        
        statmap_log( LOG_NOTICE, "Number of threads is being set to %i",  args.num_threads );
    }

    /* Set the maximum allowed reference gap between indexable subtemplates */
    if( args.assay_type == RNA_SEQ )
    {
        /* For gapped assays, set this to the default in config.h */
        args.max_reference_insert_len = REFERENCE_INSERT_LENGTH_MAX;
    } else {
        /* For ungapped assays, this should always be zero. */
        args.max_reference_insert_len = 0;
    }

    /* Set the global variables */
    num_threads = args.num_threads;
    max_reference_insert_len = args.max_reference_insert_len;
    softclip_len = args.softclip_len;

    /* Additionally soft clip by MAX_NUM_UNTEMPLATED_GS so we can go back
     * and handle them when we do the assay specific corrections */
    if( args.assay_type == CAGE )
    {
        softclip_len += MAX_NUM_UNTEMPLATED_GS;
    }

    return args;
}

/******** Read/write configuration to file ********/
void
fprintf_name_or_null( FILE* arg_fp, const char* header, char* value )
{
    fprintf( arg_fp, "%s:\t", header );
    if( NULL == value )
    {
        fprintf( arg_fp, "NULL\n" );
    } else {
        fprintf( arg_fp, "%s\n", value );
    }

    return;
}

void
write_config_file_to_stream( FILE* arg_fp, struct args_t* args  )
{
    fprintf_name_or_null( 
            arg_fp, "genome_fname", args->genome_fname );
    fprintf_name_or_null( 
            arg_fp, "genome_index_fname", args->genome_index_fname );

    fprintf_name_or_null( 
            arg_fp, "unpaired_reads_fnames", args->unpaired_reads_fnames );
    fprintf_name_or_null( 
            arg_fp, "pair1_reads_fnames", args->pair1_reads_fnames );
    fprintf_name_or_null( 
            arg_fp, "pair2_reads_fnames", args->pair2_reads_fnames );
    // struct rawread_db_t* rdb;

    fprintf_name_or_null( 
            arg_fp, "unpaired_NC_reads_fnames", args->unpaired_NC_reads_fnames );
    fprintf_name_or_null( 
            arg_fp, "pair1_NC_reads_fnames", args->pair1_NC_reads_fnames );
    fprintf_name_or_null( 
            arg_fp, "pair2_NC_reads_fnames", args->pair2_NC_reads_fnames );
    // struct rawread_db_t* NC_rdb;

    fprintf_name_or_null( 
            arg_fp, "frag_len_fname", args->frag_len_fname );
    // FILE* frag_len_fp;

    fprintf_name_or_null( 
            arg_fp, "output_directory", args->output_directory );

    fprintf_name_or_null( 
        arg_fp, "sam_output_fname", args->sam_output_fname );
    
    fprintf( arg_fp, "mapping_metaparameter:\t%.4f\n", args->mapping_metaparameter );

    fprintf( arg_fp, "num_starting_locations:\t%i\n", args->num_starting_locations );
    
    fprintf( arg_fp, "num_threads:\t%i\n", args->num_threads );

    fprintf( arg_fp, "error_model_type:\t%i\n", args->error_model_type );

    fprintf( arg_fp, "input_file_type:\t%i\n", args->input_file_type );
    fprintf( arg_fp, "assay_type:\t%i\n", args->assay_type );

    fprintf( arg_fp, "max_reference_insert_len:\t%i\n",
             args->max_reference_insert_len );

    return;
}

void
fscanf_name_or_null( FILE* arg_fp, const char* header, char** value ) 
{
    // first, build the formatting strinf
    char format_string[200];
    sprintf( format_string, "%s:\t%%s\n", header);

    // allocate space to store the read value
    char tmp_value[500];

    // scan for the value
    fscanf( arg_fp, format_string, &tmp_value );

    // if the value is equal to NULL, then set *value= NULL
    // this indicates there was no argument in the initial
    // argument list
    if( 0 == strcmp( tmp_value, "NULL" ) )
    {
        *value = NULL;
        return;
    } 
    // if it's not NULL, then allocate space for the new string,
    // and copy the option over
    else {
        *value = calloc(1, strlen(tmp_value)+1);
        strncpy( *value, tmp_value, strlen(tmp_value) );
        return;
    }
}

void
read_config_file_fname_from_disk( char* fname, struct args_t** args  )
{
    *args = malloc( sizeof(struct args_t)  );
    FILE* arg_fp = fopen( fname, "r" );

    if( arg_fp == NULL ) {
        statmap_log( LOG_FATAL, "Could not load config file '%s'",  fname );
    }

    fscanf_name_or_null( 
            arg_fp, "genome_fname", &((*args)->genome_fname) );
    fscanf_name_or_null( 
            arg_fp, "genome_index_fname", &((*args)->genome_index_fname) );

    fscanf_name_or_null( 
            arg_fp, "unpaired_reads_fnames", &((*args)->unpaired_reads_fnames) );
    fscanf_name_or_null( 
            arg_fp, "pair1_reads_fnames", &((*args)->pair1_reads_fnames) );
    fscanf_name_or_null( 
            arg_fp, "pair2_reads_fnames", &((*args)->pair2_reads_fnames) );
    // struct rawread_db_t* rdb;

    fscanf_name_or_null( 
            arg_fp, "unpaired_NC_reads_fnames", &((*args)->unpaired_NC_reads_fnames) );
    fscanf_name_or_null( 
            arg_fp, "pair1_NC_reads_fnames", &((*args)->pair1_NC_reads_fnames) );
    fscanf_name_or_null( 
            arg_fp, "pair2_NC_reads_fnames", &((*args)->pair2_NC_reads_fnames) );
    // struct rawread_db_t* NC_rdb;

    fscanf_name_or_null( 
            arg_fp, "frag_len_fname", &((*args)->frag_len_fname) );
    // FILE* frag_len_fp;

    fscanf_name_or_null( 
            arg_fp, "output_directory", &((*args)->output_directory) );

    fscanf_name_or_null( 
        arg_fp, "sam_output_fname", &((*args)->sam_output_fname) );
    
    fscanf( arg_fp, "mapping_metaparameter:\t%f\n",
            &((*args)->mapping_metaparameter) );

    fscanf( arg_fp, "num_starting_locations:\t%i\n", 
            &((*args)->num_starting_locations) );
    
    fscanf( arg_fp, "num_threads:\t%i\n", 
            &((*args)->num_threads) );

    fscanf( arg_fp, "error_model_type:\t%u\n",
            &((*args)->error_model_type) );

    fscanf( arg_fp, "input_file_type:\t%u\n", 
            &( (*args)->input_file_type) );
    fscanf( arg_fp, "assay_type:\t%u\n", 
            &( (*args)->assay_type) );

    fscanf( arg_fp, "max_reference_insert_len:\t%i\n",
            &( (*args)->max_reference_insert_len ) );

    fclose( arg_fp  );

    return;
}

void
read_config_file_from_disk( struct args_t** args  )
{
    read_config_file_fname_from_disk( "config.dat", args );
}
