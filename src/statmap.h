/* Copyright (c) 2009-2010 Nathan Boley */

extern int num_threads;
extern int min_num_hq_bps;

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

