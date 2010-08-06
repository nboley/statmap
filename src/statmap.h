/* Copyright (c) 2009-2010 Nathan Boley */

#ifndef STATMAP_H
#define STATMAP_H

#include "config.h"

struct genome_data;

extern int num_threads;
extern int min_num_hq_bps;

/*
 * The assay that generated these reads.
 * 
 * This only matters if we want to do 
 * iterative mapping, and for the wiggle
 * output type.
 *
 */
enum assay_type_t {
    // UNKNOWN = 0,
    CAGE = 1,
    CHIP_SEQ = 2
};

extern int num_trace_tracks;
extern char** trace_track_names;


typedef enum {
    // We have not specified the file type
    // UNKNOWN = 0,
    // Fastq file. Quality scores are PHRED 
    // ( offset 40, log10 prob )
    SANGER_FQ = 1,
    // Illumina version 1.3
    // Fastq file. Quality scores are PHRED 
    // ( offset 64, log10 prob )
    ILLUMINA_v13_FQ = 2,
    // Solexa version 1.4+
    // Fastq file. Quality scores are PHRED 
    // ( offset 64, log10 prob )
    SOLEXA_v14_FQ = 2,
    // Solexa version 1.5+
    // Fastq file. Quality scores are PHRED 
    // ( offset 64, log10 prob )
    // In addition, scores of 0 and 1 are never used
    // Furthermore, a score of 2 is special: 
    /* 
     *  "A quality score of 2, encoded as a "B", is used as a special indicator. A
     *   quality score of 2 does not imply a specific error rate, but rather implies that
     *   the marked region of the read should not be used for downstream analysis."
     *   http://docs.google.com/fileview?id=0B-lLYVUOliJFYjlkNjAwZjgtNDg4ZC00MTIyLTljNjgtMmUzN2M0NTUyNDE3&hl=en&pli=1
     */
    ILLUMINA_v15_FQ = 2,
    // Solexa version 1.0
    // Fastq file. Quality scores are 
    // ( offset 64, log10 log odds )
    SOLEXA_LOG_ODDS_FQ = 3,
    // Test suite error format
    // All scores are perfect h's, 
    // assume it's the test code, and the 
    // scores are ILLUMINA_v13_FQ
    TEST_SUITE_FORMAT = 2,
    // Marks SOLEXA ( from his ChIP-seq )
    // Maximum is ASCII V
    MARKS_SOLEXA = 4
} input_file_type_t;

/* fwd declaration for the rawread db type */
struct rawread_db_t;

/* Store parsed options */
typedef struct {
    struct genome_data* genome;
    char* genome_fname;
    char* genome_index_fname;
    
    char* unpaired_reads_fnames;
    char* pair1_reads_fnames;
    char* pair2_reads_fnames;
    struct rawread_db_t* rdb;

    char* unpaired_NC_reads_fnames;
    char* pair1_NC_reads_fnames;
    char* pair2_NC_reads_fnames;
    struct rawread_db_t* NC_rdb;

    char* snpcov_fname;
    FILE* snpcov_fp;

    char* frag_len_fname;
    FILE* frag_len_fp;

    char* output_directory;

    char* sam_output_fname;

    char* log_fname;
    FILE* log_fp;

    float min_match_penalty;
    float max_penalty_spread;
    int min_num_hq_bps;

    int num_starting_locations;

    int indexed_seq_len;
    int num_threads;
        
    input_file_type_t input_file_type;
    enum assay_type_t assay_type; 
} args_t;

void usage();

/*
 * Try and determine the file type. 
 *
 * In particular, we try and determine what the sequence mutation
 * string types are.
 *
 * The method is to scan the first 10000 reads and record the 
 * max and min untranslated scores. 
 *
 */

input_file_type_t
guess_input_file_type( args_t* args );

/*
 * Guess the optimal indexed sequence length.
 *
 * To do this, we open up the first read file and then scan for a read. The 
 * seq length of the first read is what we set the index length to.
 *
 */
int
guess_optimal_indexed_seq_len( args_t* args);

args_t
parse_arguments( int argc, char** argv );

struct mapped_reads_db;

void
map_marginal( args_t* args, 
              struct genome_data* genome,
              struct rawread_db_t* rdb,
              struct mapped_reads_db** mpd_rds_db,
              enum bool is_nc );

void
build_fl_dist( args_t* args, struct mapped_reads_db* mpd_rds_db );

void
iterative_mapping( args_t* args, 
                   struct genome_data* genome,
                   struct mapped_reads_db* mpd_rds_db );

int 
main( int argc, char** argv );


#endif /* STATMAP_H */
