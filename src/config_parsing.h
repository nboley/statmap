/* Copyright (c) 2009-2010, Nathan Boley */

#ifndef CONFIG_PARSING_H
#define CONFIG_PARSING_H

#include <stdio.h>
#include "config.h"

enum input_file_type_t {
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
     *   quality score of 2 does not imply a specific error rate, but rather implies 
     *   that the marked region of the read should not be used for downstream analysis."
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
};

/* Store parsed options */
struct args_t {
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
        
    enum input_file_type_t input_file_type;
    enum assay_type_t assay_type; 
};

void
write_config_file_to_disk( struct args_t* args  );

/* this assumes that we have moved intot he output directory */
void
read_config_file_from_disk( struct args_t** args  );

#endif // CONFIG_PARSING_H
