/* Copyright (c) 2009-2010, Nathan Boley */

#include <pthread.h>
#include <stdio.h>

#include "genome.h"
#include "rawread.h"
#include "db_interface.h" 

struct single_map_thread_data {
    int thread_id;
    struct genome_data* genome;
    FILE* log_fp;
    pthread_mutex_t* log_fp_mutex;

    struct rawread_db_t* rdb;
    
    unsigned int* mapped_cnt;
    pthread_mutex_t* mapped_cnt_mutex;
    
    candidate_mappings_db* mappings_db;
    pthread_mutex_t* mappings_db_mutex;
    float min_match_penalty;
    float max_penalty_spread;
    int max_subseq_len;
};


void*
find_candidate_mappings( void* params );


void
find_all_candidate_mappings( struct genome_data* genome,
                             FILE* log_fp,
                             struct rawread_db_t* rdb,

                             candidate_mappings_db* mappings_db,
                             float min_match_penalty,
                             float max_penalty_spread,
                             float max_seq_length

    );
