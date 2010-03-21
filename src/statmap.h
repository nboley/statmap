/* Copyright (c) 2009-2010 Nathan Boley */

#ifndef STATMAP_H
#define STATMAP_H

#include "genome.h"

extern int num_threads;
extern int min_num_hq_bps;

/*
 * The assay that generated these reads.
 * 
 * This only matters if we want to do 
 * iterative mapping.
 *
 */
enum assay_type_t {
    // UNKNOWN = 0,
    CAGE = 1,
    CHIP_SEQ = 2
};


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
    SOLEXA_v14p_FQ = 2,
    // Solexa version 1.0
    // Fastq file. Quality scores are 
    // ( offset 64, log10 log odds )
    SOLEXA_LOG_ODDS_FQ = 3,
    // Test suite error format
    // All scores are perfect h's, 
    // assume it's the test code, and the 
    // scores are ILLUMINA_v13_FQ
    TEST_SUITE_FORMAT = 2

} input_file_type_t;

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

    int num_threads;
    int min_num_hq_bps;
    
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

/* 
   
TODO - MOVE THESE
   
void
join_all_candidate_mappings( 
    candidate_mappings_db* cand_mappings_db,
    mapped_reads_db* mpd_rds_db );
   
void
write_mapped_reads_to_sam( rawread_db_t* rdb,
                           mapped_reads_db* mappings_db,
                           genome_data* genome,
                           FILE* sam_ofp );
// END MOVE THESE
   
*/ 

void
map_marginal( args_t* args, genome_data* genome );

int 
main( int argc, char** argv );


#endif /* STATMAP_H */
