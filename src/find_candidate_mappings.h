/* Copyright (c) 2009-2010, Nathan Boley */

#include <pthread.h>
#include <stdio.h>

#include "genome.h"
#include "rawread.h"
#include "candidate_mapping.h"
#include "mapped_location.h"

void
search_index( struct index_t* index, 
              
              float min_match_penalty,
              float max_penalty_spread,
              mapped_locations** results,

              struct rawread* r,
              float* bp_mut_rates,

              float* lookuptable_position,
              float* inverse_lookuptable_position,
              float* reverse_lookuptable_position,
              float* reverse_inverse_lookuptable_position
    );


struct single_map_thread_data {
    int thread_id;
    struct genome_data* genome;
    FILE* log_fp;
    pthread_mutex_t* log_fp_mutex;

    struct rawread_db_t* rdb;
    
    unsigned int* mapped_cnt;
    pthread_mutex_t* mapped_cnt_mutex;
    readkey_t max_readkey;
    
    candidate_mappings_db* mappings_db;
    pthread_mutex_t* mappings_db_mutex;
    float min_match_penalty;
    float max_penalty_spread;
    int max_subseq_len;

    /* Pointer to error data struct shared by all threads */
    struct error_data_t* global_error_data;
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
