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

              struct penalty_array_t* pa,

              bool only_find_unique_mappers
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

    /* global_error_data contains the weighted average of error data so far
     * and is used to calculate the lookup tables */
    struct error_data_t* global_error_data;
    /* scratch_error_data contains the error data from the current set of
     * running threads, and is synchronized with global_error_data every
     * READS_STAT_UPDATE_STEP_SIZE reads */
    struct error_data_t* scratch_error_data;

    /* fn ptr to build penalty array for a read */
    void (*build_penalty_array_from_rawread)(
            struct rawread*,
            struct error_data_t*,
            struct penalty_array_t*
        );

    bool only_find_unique_mappers;
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
                             float max_seq_length,

                             enum SEARCH_TYPE search_type
    );

