/* Copyright (c) 2009-2010, Nathan Boley */

#ifndef CONFIG_PARSING_H
#define CONFIG_PARSING_H

#include <stdio.h>
#include "config.h"

/* Store parsed command line options */
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

    char* frag_len_fname;
    FILE* frag_len_fp;

    char* output_directory;

    char* sam_output_fname;
    
    float min_match_penalty;
    float max_penalty_spread;
    int min_num_hq_bps;

    int num_starting_locations;
    
    int num_threads;

    enum SEARCH_TYPE search_type;
        
    enum input_file_type_t input_file_type;
    enum assay_type_t assay_type; 
};

struct args_t
parse_arguments( int argc, char** argv );

void
write_config_file_to_stream( FILE* arg_fp, struct args_t* args  );

/* this assumes that we have moved intot he output directory */
void
read_config_file_fname_from_disk( char* fname, struct args_t** args  );

void
read_config_file_from_disk( struct args_t** args  );

#endif // CONFIG_PARSING_H
